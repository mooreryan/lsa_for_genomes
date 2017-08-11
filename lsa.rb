#!/usr/bin/env ruby

PLOT_FUNCTION = %Q{
## This cutoff index is 1-based
plot.colored.by.inflection.point <- function(dat, cutoff, xlab="Rank", ylab="Weight")
{
    par(mfrow=c(1,1), lwd=2, cex.axis=0.8, cex.lab=1.2)

    col1 <- rgb(1, 0, 0, 0.5)
    col2 <- rgb(0, 0, 0, 0.5)

    plot(dat,
         xlab="", ylab="",
         bty="n",
         axes=F,
         type="n")

    grid(lwd=1)
    box()
    axis(1)
    axis(2, las=1)
    title(xlab=xlab, line=2.5)
    title(ylab=ylab, line=2.75)

    points(x=1:cutoff,
           y=dat[1:cutoff],
           col=col1,
           pch=16,
           cex=0.8)
    points(x=cutoff+1:length(dat),
           y=dat[cutoff+1:length(dat)],
           col=col2,
           pch=16,
           cex=0.8)
}
}

Signal.trap("PIPE", "EXIT")

require "tempfile"

require "aai"
require "abort_if"
require "fileutils"
require "matrix"
require "set"
require "trollop"

require "lsa"

include AbortIf
Process.extend Aai::CoreExtensions::Process

module Aai
  extend Aai
  extend Aai::Utils
end

module Lsa
  PIPELINE_VERSION = "0.5.0"
  COPYRIGHT = "2017 Ryan Moore"
  CONTACT   = "moorer@udel.edu"
  WEBSITE   = "https://github.com/mooreryan/lsa_for_genomes"
  LICENSE   = "MIT"

  VERSION_BANNER = "  # Pipeline version: #{PIPELINE_VERSION}
  # Lib version: #{VERSION}
  # Copyright #{COPYRIGHT}
  # Contact: #{CONTACT}
  # Website: #{WEBSITE}
  # License: #{LICENSE}"

  extend Lsa
end

def abort_unless_command exe
  abort_unless_file_exists exe

  abort_if File.directory?(exe),
           "#{exe} is a directory, should be a command"

  abort_unless File.executable?(exe),
               "File #{exe} is not executable"
end

def all_files_exist? fnames
  fnames.all? { |fname| File.exists? fname }
end

def clean_str str
  str.strip.gsub(/[^\p{Alnum}_]+/, "_").gsub(/_+/, "_")
end

def make_dirs *dirs
  dirs.each do |dir|
    FileUtils.mkdir_p dir
  end
end

def parse_index_mapping_file fname
  idx2data = {}

  File.open(fname, "rt").each_line do |line|
    idx, data = line.chomp.split "\t"

    idx2data[idx.to_i] = data
  end

  idx2data
end

def max_topics num_terms, num_docs
  [num_terms, num_docs].min
end

def run_redsvd redsvd,
               td_matrix_outfname,
               redsvd_outf_basename,
               num_topics,
               log_out_fname,
               log_err_fname

  cmd = "#{redsvd} " +
        "--format sparse " +
        "--method SVD " +
        "--input #{td_matrix_outfname} " +
        "--output #{redsvd_outf_basename} " +
        "--rank #{num_topics} " +
        "1>> #{log_out_fname} 2>> #{log_err_fname}"

  Process.run_and_time_it! "Running the LSA transform", cmd
end

def run_td_matrix_counts td_matrix,
                         cluster_outfname,
                         outdir,
                         log_out_fname,
                         log_err_fname,
                         td_matrix_outfiles

  if all_files_exist? td_matrix_outfiles
    AbortIf.logger.info { "Term-doc matrix already created, skipping" }
  else
    cmd = "#{td_matrix} " +
          "#{cluster_outfname} " +
          "#{outdir} " +
          "1>> #{log_out_fname} 2>> #{log_err_fname}"
    Process.run_and_time_it! "Building term-doc matrix", cmd
  end
end

def vector_projection a, b
  (a.dot(b) / b.dot(b).to_f) * b
end

def vector_rejection a, b
  a - vector_projection(a, b)
end

def index_of_max ary
  ary.each_with_index.max.last
end

# Returns the index of the inflection point
def inflection_point dat
  num_points = dat.length
  first_point = Vector[0, dat.first]
  last_point = Vector[num_points-1, dat.last]

  b = last_point - first_point

  dists = (1..num_points-2).map do |idx|
    this_point = Vector[idx, dat[idx]]
    a = this_point - first_point
    vector_rejection(a, b).norm
  end

  index_of_max(dists) + 1 # add one for the 1.. range
end

THIS_DIR = File.absolute_path(File.join File.dirname __FILE__)

opts = Trollop.options do
  version Lsa::VERSION_BANNER

  banner <<-EOS

#{Lsa::VERSION_BANNER}

  Latent semantic analysis pipeline for genomes and metagenomes.

  --num-topics is the LIMIT on topics. You may see fewer depending on
    the data. If you pass --num-topics 0, all topics will be used.

  --percent-of-terms-per-topic has an automatic mode, pass 0

  Options:
  EOS

  opt(:binary_dir,
      "Folder with the LSA scripts and binaries",
      default: File.join(THIS_DIR, "bin"))
  opt(:mmseqs,
      "Location of the mmseqs binary",
      default: File.join(Dir.home, "bin", "mmseqs"))

  opt(:infiles,
      "Files with ORFs",
      type: :strings)
  opt(:outdir,
      "Output directory",
      type: :string,
      default: File.join(THIS_DIR, "output"))
  opt(:mapping,
      "Mapping file (optional)",
      type: :string)

  opt(:cpus,
      "Number of CPUs to use",
      default: 4)

  opt(:num_topics,
      "The maximum number of topics to calculate for LSA" +
      " (Use 0 for maximum number of topics.)",
      default: 0)
  opt(:percent_of_terms_per_topic,
      "What percentage of top terms per topic do you want to look at?" +
      " (Use 0 for automatic)",
      default: 0)
end

abort_unless opts[:percent_of_terms_per_topic] >= 0 &&
             opts[:percent_of_terms_per_topic] <= 100,
             "--percent-of-terms-per-topic must be >= 0 and <= 100"

abort_unless Dir.exist?(opts[:binary_dir]),
             "The directory specified by --bin-dir doesn't exist"

######################################################################
# check commands
################


prep_seq_files = File.join opts[:binary_dir], "prep_seq_files.rb"
cluster = File.join opts[:binary_dir], "cluster.rb"
td_matrix = File.join opts[:binary_dir], "td_matrix"
redsvd = File.join opts[:binary_dir], "redsvd"
mmseqs = opts[:mmseqs]
make_color_maps = File.join opts[:binary_dir], "make_color_maps.rb"

abort_unless_command prep_seq_files
abort_unless_command cluster
abort_unless_command td_matrix
abort_unless_command redsvd
abort_unless_command mmseqs
abort_unless_command make_color_maps

################
# check commands
######################################################################

######################################################################
# outfiles
##########

seq_dir = File.join opts[:outdir], "prepped_seqs"
cluster_dir = File.join opts[:outdir], "clustering"
log_dir = File.join opts[:outdir], "logs"

FileUtils.mkdir_p seq_dir
FileUtils.mkdir_p cluster_dir
FileUtils.mkdir_p log_dir

prepped_seq_files_basename = "seqs"
prepped_seq_files = File.join seq_dir, prepped_seq_files_basename + ".faa"

mmseqs_base_with_dir = File.join cluster_dir, prepped_seq_files_basename

mmseqs_final_outfname = "#{mmseqs_base_with_dir}.clu.tsv.sorted"
mmseqs_outfiles = [
  "#{mmseqs_base_with_dir}.DB",
  "#{mmseqs_base_with_dir}.DB.index",
  "#{mmseqs_base_with_dir}.DB.lookup",
  "#{mmseqs_base_with_dir}.DB_h",
  "#{mmseqs_base_with_dir}.DB_h.index",
  "#{mmseqs_base_with_dir}.clu",
  "#{mmseqs_base_with_dir}.clu.index",
  "#{mmseqs_base_with_dir}.clu.tsv",
  mmseqs_final_outfname,
]


log_out_fname = File.join log_dir, "lsa_out.txt"
log_err_fname = File.join log_dir, "lsa_err.txt"

##########
# outfiles
######################################################################


######################################################################
# check files, prep dirs
########################

abort_if opts[:infiles].nil? || opts[:infiles].empty?,
         "No infiles given. Try #{__FILE__} --help for help."

# abort_if_file_exists opts[:outdir] unless opts[:force]

FileUtils.mkdir_p opts[:outdir]

abort_if opts[:mapping] && !File.exist?(opts[:mapping]),
         "--mapping #{opts[:mapping]} was passed, but " +
         "#{opts[:mapping]} does not exist."

########################
# check files, prep dirs
######################################################################

######################################################################
# run the pipeline
##################


# Prep seq files

if File.exists? prepped_seq_files
  AbortIf.logger.info { "Seqs already prepped, skipping" }
else
  infiles = opts[:infiles].join " "
  cmd = "#{prep_seq_files} #{infiles} 1> #{prepped_seq_files} 2> #{log_err_fname}"
  Process.run_and_time_it! "Prepping ORFs", cmd
  abort_unless_file_exists prepped_seq_files
end

##########
##########
##########
# Cluster

if all_files_exist? mmseqs_outfiles
  AbortIf.logger.info { "Clustering already done, skipping" }
else
  cmd = "#{cluster} #{mmseqs} #{opts[:cpus]} #{prepped_seq_files} #{cluster_dir} #{prepped_seq_files_basename} 1>> #{log_out_fname} 2>> #{log_err_fname}"
  Process.run_and_time_it! "Clustering ORFs", cmd
end

##########
##########
##########

label2outf = nil
doc2new_doc = nil
new_cluster_outfnames = []
color_map_tmp_dir = nil
if opts[:mapping]
  # Read the mapping file
  AbortIf.logger.info { "Parsing mapping file" }
  label2outf, doc2new_doc = Lsa::parse_mapping_file opts[:mapping],
                                                    mmseqs_final_outfname

  # Make the new cluster files
  AbortIf.logger.info { "Making new cluster files" }
  new_cluster_outfnames =
    Lsa::make_new_cluster_files mmseqs_final_outfname,
                                label2outf,
                                doc2new_doc

  # Make color maps for iroki
  color_map_dir = File.join opts[:outdir], "color_maps"
  cmd = [make_color_maps, opts[:mapping], color_map_dir].join " "
  Process.run_and_time_it! "Making color maps", cmd
end

##########
##########
##########

# From here on, run everything in a big loop for all the cluster files
[mmseqs_final_outfname, new_cluster_outfnames].
  flatten.each_with_index do |cluster_outfname, metadata_group_idx|

  AbortIf.logger.info { "Running big loop on #{cluster_outfname}" }

  if metadata_group_idx.zero?
    metadata_group_label = "original"
  else
    metadata_group_label = label2outf.to_a[metadata_group_idx-1].first
  end

  metadata_group_dir =
    File.join opts[:outdir], "metadata_groups", metadata_group_label
  td_matrix_dir =
    File.join metadata_group_dir, "td_matrix"
  redsvd_dir =
    File.join metadata_group_dir, "redsvd"
  rscript_dir =
    File.join metadata_group_dir, "r"
  trees_dir =
    File.join metadata_group_dir, "trees"
  top_terms_by_topic_dir =
    File.join metadata_group_dir, "top_terms_by_topic"
  top_terms_by_doc_dir =
    File.join metadata_group_dir, "top_terms_by_doc"

  make_dirs metadata_group_dir,
            td_matrix_dir,
            redsvd_dir,
            rscript_dir,
            trees_dir,
            top_terms_by_topic_dir,
            top_terms_by_doc_dir

  td_matrix_outfname =
    File.join td_matrix_dir, "td_matrix.txt"
  idx_to_doc_outfname =
    File.join td_matrix_dir, "td_matrix.idx_to_doc_map.txt"
  idx_to_term_outfname =
    File.join td_matrix_dir, "td_matrix.idx_to_term_map.txt"

  td_matrix_outfiles = [
    td_matrix_outfname,
    idx_to_doc_outfname,
    idx_to_term_outfname,
  ]

  redsvd_outf_base =
    File.join redsvd_dir, "svd"

  make_trees_r_script_fname =
    File.join rscript_dir, "make_tree.r"

  plot_by_topic_term_weights_r_script_fname =
    File.join rscript_dir, "by_topic_term_weights.r"


  newick_docs_fname =
    File.join trees_dir, "doc_dist_tree.newick.txt"

  svd_U_fname            = "#{redsvd_outf_base}.U"
  svd_S_fname            = "#{redsvd_outf_base}.S"
  svd_V_fname            = "#{redsvd_outf_base}.V"
  svd_VS_fname           = "#{redsvd_outf_base}.VS"
  svd_US_fname           = "#{redsvd_outf_base}.US"
  svd_VS_dis_fname       = "#{redsvd_outf_base}.VS.dis"
  svd_VS_dis_fname_for_r = "#{redsvd_outf_base}.VS.dis.for_r"
  svd_VS_US_dis_fname    = "#{redsvd_outf_base}.VS_to_US.dis"

  top_terms_by_topic =
    File.join top_terms_by_topic_dir, "top_terms_by_topic.txt"
  top_terms_by_topic_plot_basename =
    File.join top_terms_by_topic_dir, "top_terms_plot"

  top_terms_by_doc =
    File.join top_terms_by_doc_dir, "top_terms_by_doc.txt"
  top_terms_by_doc_plot_basename =
    File.join top_terms_by_doc_dir, "top_terms_plot"

  run_td_matrix_counts td_matrix,
                       cluster_outfname,
                       File.dirname(td_matrix_outfname),
                       log_out_fname,
                       log_err_fname,
                       td_matrix_outfiles

  # If we are in the original metadata group, move the mapping files
  # into the trees dir
  if opts[:mapping] && metadata_group_label == "original"
    FileUtils.mv color_map_dir, trees_dir
  end

  # Parse the idx to name map files
  idx2doc  = parse_index_mapping_file idx_to_doc_outfname
  idx2term = parse_index_mapping_file idx_to_term_outfname

  num_docs  = idx2doc.count
  num_terms = idx2term.count

  num_terms_to_keep = nil
  sorted_weights_with_index = nil

  # If the user didn't specify number of topics to use, select the max
  # possible topics. I.e., the smaller of number of terms and number
  # of docs.
  if opts[:num_topics].zero?
    opts[:num_topics] = max_topics num_terms, num_docs
  end

  # LSA
  run_redsvd redsvd,
             td_matrix_outfname,
             redsvd_outf_base,
             opts[:num_topics],
             log_out_fname,
             log_err_fname

  # Make the R script for the tree

  AbortIf.logger.info do
    "Finding top terms for each topic (weight listed is the raw " +
      "weight not the projection)"
  end

  # Read the matrix
  topic2top_terms = {}
  topics = File.open(svd_US_fname, "rt").
           read.
           chomp.
           split("\n").
           map { |line| line.split(" ").map(&:to_f) }.
           transpose

  # For the R graphs
  num_terms_to_keep_per_topic = []
  topic_plot_fnames = []

  # Write top terms
  File.open(top_terms_by_topic, "w") do |f|
    topics.each_with_index do |topic, topic_idx|
      sorted_weights_with_index = topic.map.
                                  with_index do |weight, term_idx|
        abort_unless idx2term.has_key?(term_idx),
                     "#{term_idx} missing from idx2term"

        term_name = idx2term[term_idx]

        [term_name, weight, weight.abs]
      end.sort_by do |term, weight, abs_weight|
        abs_weight
      end.reverse

      if opts[:percent_of_terms_per_topic].zero?
        sorted_abs_weights = sorted_weights_with_index.map(&:last)
        num_terms_to_keep = inflection_point(sorted_abs_weights) + 1
      else
        num_terms_to_keep = (num_terms * (opts[:percent_of_terms_per_topic] / 100.0)).round
      end
      num_terms_to_keep_per_topic << num_terms_to_keep
      topic_plot_fnames <<
        "#{top_terms_by_topic_plot_basename}." +
        "topic_#{topic_idx}.pdf"

      AbortIf.logger.info do
        "Keeping #{num_terms_to_keep} of #{num_terms} terms for " +
          "topic #{topic_idx}"
      end

      top_terms = sorted_weights_with_index.
                  take(num_terms_to_keep)

      # Keep the term names in the ht
      topic2top_terms[topic_idx] = Set.new(top_terms.map(&:first))

      top_terms.each do |term, weight, abs_weight|
        f.puts [topic_idx, term, weight].join " "
      end
    end
  end

  # Plot the term weights per topic plots
  rscript_str = %Q{
#{PLOT_FUNCTION}

dat <- read.table("#{svd_US_fname}", sep=" ")

cutoffs <- c(#{num_terms_to_keep_per_topic.join(", ")})
plot_fnames <- c(#{topic_plot_fnames.map{|str| %Q["#{str}"]}.join(", ")})

for (cidx in 1:ncol(dat)) {
    pdf(plot_fnames[cidx], width=8, height=5)
    plot.colored.by.inflection.point(sort(abs(dat[, cidx]), decreasing=T), cutoffs[cidx])
    invisible(dev.off())
}

}

  File.open(plot_by_topic_term_weights_r_script_fname, "w") do |f|
    f.puts rscript_str
  end

  cmd = "Rscript #{plot_by_topic_term_weights_r_script_fname}"
  Process.run_and_time_it! "Ploting by topic term weights", cmd

  ##########
  ##########
  ##########
  ##########

  doc_names = []
  AbortIf.logger.info { "Finding top terms for each document" }
  term2cluster = {}
  # This is for coloring the points for each of the graphs in R
  num_terms_to_keep_per_document = []

  # Inflection point is not great for term to doc weights as they
  # generally cluster in groups
  File.open(top_terms_by_doc, "w") do |top_terms_f|
    top_terms_f.puts %w[doc_name doc_idx term_name term_idx cluster dist].join " "
    # Get terms closest to each doc. Each line is the dists to each
    # term for that doc.
    File.open(svd_VS_US_dis_fname, "rt").each_line.with_index do |line, doc_idx|
      # Need some output files
      abort_unless idx2doc.has_key?(doc_idx),
                   "Doc index #{doc_idx} is missing from the idx2doc hash table"
      doc_name = idx2doc[doc_idx]
      doc_names << doc_name

      doc_dir = File.join top_terms_by_doc_dir, doc_name
      FileUtils.mkdir_p doc_dir

      dists_fname =
        File.join doc_dir,
                  "term_doc_dists.txt"
      clusters_fname =
        File.join doc_dir,
                  "term_doc_dists_clusters.txt"
      kmeans_r_fname =
        File.join rscript_dir,
                  "term_doc_dists_kmeans.r"
      centroids_plot_fname =
        File.join doc_dir,
                  "term_doc_dists_cluster_centroids_plot.pdf"
      dists_plot_fname =
        File.join doc_dir,
                  "term_doc_dists_plot.pdf"

      *dists = line.chomp.split " "

      File.open(dists_fname, "w") do |f|
        dists.each { |dist| f.puts dist }
      end

      abort_unless dists.count == idx2term.count,
                   "Number of terms from #{svd_VS_US_dis_fname} does not match number of terms in idx2term hash table."

      dists_with_idx =
        dists.map.with_index { |dist, term_idx| [term_idx, dist] }

      dist_counts =
        dists.group_by(&:itself).map { |dist, ary| [dist, ary.count] }.sort_by { |dist, count| count}.reverse

      # No + 1 because we want one less than the inflection point
      num_centroids = inflection_point(dist_counts.map(&:last))
      num_centroids = 2 if num_centroids == 1
      centroids = dist_counts.take(num_centroids).map(&:first)
      AbortIf.logger.info { "#{doc_name} will have #{centroids.count}" }

      # START HERE need to make this R code good, and change the above file names
      rscript_str = %Q{
#{PLOT_FUNCTION}

num.centroids <- #{num_centroids}
dat <- read.table("#{dists_fname}", col.names=c("dist"))
dat.sorted <- dat[with(dat, order(dist)),]
dist.counts <- c(#{dist_counts.map { |dist, count| count }.join(", ")})

## Get the clusters
k <- kmeans(dat.sorted, c(#{centroids.join(", ")}))
write.table(k$cluster, "#{clusters_fname}", row.names=F, col.names=F, quote=F)

## Plot dist counts to show how number of centroids was determined
## TODO should make this more tolerant like call a mean dist equal if
## within a certain tolerance?
pdf("#{centroids_plot_fname}", width=8, height=5)
plot.colored.by.inflection.point(dist.counts, num.centroids, ylab="Distance")
invisible(dev.off())

## Plot how the data fits into the clusters
pdf("#{dists_plot_fname}", width=8, height=5)
par(mfrow=c(1,1), lwd=2, cex.axis=0.8, cex.lab=1.2)
plot(dat.sorted, xlab="", ylab="", bty="n", axes=F, type="n")
grid(lwd=1)
box()
axis(1)
axis(2, las=1)
title(xlab="Rank", line=2.5)
title(ylab="Distance", line=2.75)
points(dat.sorted, col=rainbow(num.centroids)[k$cluster], pch=16, cex=0.8)
invisible(dev.off())
}

      File.open(kmeans_r_fname, "w") do |f|
        f.puts rscript_str
      end

      Process.run_and_time_it! "Running kmeans", "Rscript #{kmeans_r_fname}"

      clusters = File.open(clusters_fname, "rt").read.chomp.split("\n")

      dists_w_cluster =
        dists.map.with_index { |dist, idx| [clusters[idx], dist, idx] }

      dists_w_cluster.each do |cluster, dist, term_idx|
        abort_unless idx2term.has_key?(term_idx),
                     "Missing term index #{term_idx} from the idx2term hash table."

        term_name = idx2term[term_idx]

        # TODO this might break if the seq headers are not unique
        # across the different input files
        #
        # Need to key the term names not on the whole header, but the
        # part after the tilde, because the part in front of the tilde
        # (the doc part) changes for the non-original metadata
        # categories
        term_name_no_tilde = term_name.split("~").last

        # The term2cluster hash table...each term will have a cluster
        # number for each doc.
        if term2cluster.has_key? term_name_no_tilde
          abort_if term2cluster[term_name_no_tilde].has_key?(doc_name),
                   "Term #{term_name_no_tilde} was repeated for doc #{doc_name} in term2cluster hash table"

          term2cluster[term_name_no_tilde][doc_name] = cluster
        else
          term2cluster[term_name_no_tilde] = { doc_name => cluster }
        end

        top_terms_f.puts [doc_name, doc_idx, term_name_no_tilde, term_idx, cluster, dist].join " "
      end
    end
  end

  abort_unless term2cluster.count == idx2term.count,
               "There were #{idx2term.count} terms in idx2term, but term2cluster only has #{term2cluster.count} terms"

  ##########
  ##########
  ##########
  ##########


  AbortIf.logger.info { "Converting VS.dis file to R format" }
  # Convert the VS.dis file to an R ready format
  nrows = 0
  mat = nil
  File.open(svd_VS_dis_fname, "rt").each_line.with_index do |line, idx|
    if idx.zero?
      nrows = line.chomp.to_i
      mat = Array.new(nrows) { Array.new(nrows, 0) }
    else
      ridx, cidx, val = line.chomp.split " "
      ridx = ridx.to_i
      cidx = cidx.to_i
      val = val.to_f

      mat[ridx][cidx] = val
      mat[cidx][ridx] = val
    end
  end

  # Write out the dist file in a format that R can read
  File.open(svd_VS_dis_fname_for_r, "w") do |f|
    abort_unless idx2doc.count == mat.count,
                 "Size of dist matrix doesn't match number of docs in idx2doc"

    # write the header line
    f.puts idx2doc.sort_by { |idx, doc| idx }.map(&:last).join " "

    mat.each do |row|
      f.puts row.join " "
    end
  end

  rscript_str = %Q{
library("ape")

proj.docs.dist <- as.dist(read.table("#{svd_VS_dis_fname_for_r}", header=T, sep=" "))

proj.docs.dist.tree <- as.phylo(hclust(proj.docs.dist, method="average"))

write.tree(proj.docs.dist.tree, file="#{newick_docs_fname}")
}

  File.open(make_trees_r_script_fname, "w") do |f|
    f.puts rscript_str
  end

  cmd = "Rscript #{make_trees_r_script_fname}"
  Process.run_and_time_it! "Generating trees", cmd

  top_terms_per_topic_fnames = {}
  topic2top_terms.each do |topic, terms|
    dir = File.join top_terms_by_topic_dir, "seqs"
    FileUtils.mkdir_p dir
    fname = File.join dir,
                      "top_terms.topic_#{topic.to_i}.fa"

    top_terms_per_topic_fnames[topic] = File.open fname, "w"
  end

  term_doc_dist_cluster_fnames = {}
  doc_names.each do |doc_name|
    abort_if term_doc_dist_cluster_fnames.has_key?(doc_name),
             "Doc #{doc_name} is repeated in term_doc_dist_cluster_fnames hash table"

    # This will need to be filled out as you read the cluster file?
    term_doc_dist_cluster_fnames[doc_name] = {}
  end

  AbortIf.logger.info { "Grepping seq ids" }
  # Grep them from the prepped file
  ParseFasta::SeqFile.open(prepped_seq_files).each_record do |rec|
    # If the sequence is in this file but not in term2cluster hash
    # table, then that sequence was NOT a cluster rep seq from MMseqs2
    # clustering

    header_no_tilde = rec.header.split("~").last
    rec.header = header_no_tilde

    if term2cluster.has_key? header_no_tilde
      term2cluster[header_no_tilde].each do |doc_name, cluster|
        abort_unless term_doc_dist_cluster_fnames.has_key?(doc_name),
                     "Doc #{doc_name} present in term2cluster[header_no_tilde] but missing from term_doc_dist_cluster_fnames hash table"

        unless term_doc_dist_cluster_fnames[doc_name].has_key?(cluster)
          doc_dir = File.join top_terms_by_doc_dir, doc_name
          abort_unless File.exist?(doc_dir),
                       "#{doc_dir} does not exist, but it should have already been created"
          seqs_dir = File.join doc_dir, "seqs"
          FileUtils.mkdir_p seqs_dir

          fname = File.join seqs_dir, "#{doc_name}_term_doc_dist_cluster_#{cluster}.fa"
          term_doc_dist_cluster_fnames[doc_name][cluster] = File.open fname, "w"
        end

        term_doc_dist_cluster_fnames[doc_name][cluster].puts rec
      end
    end

    topic2top_terms.each do |topic, headers|
      # TODO this will break if headers are not unique across original files?
      headers_no_tilde = Set.new(headers.map{|header| header.split("~").last})

      if headers_no_tilde.include? header_no_tilde
        # seqs can be top seq in multiple topics
        top_terms_per_topic_fnames[topic].puts rec
      end
    end
  end

  # Close the files
  top_terms_per_topic_fnames.each do |topic, f|
    f.close
  end

  top_terms_per_doc_fnames.each do |doc, f|
    f.close
  end
end

##################
# run the pipeline
######################################################################

AbortIf.logger.info { "Done!" }

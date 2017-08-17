#!/usr/bin/env ruby

def biplot_rscript svd_US_fname, svd_VS_fname, doc_names, pdf_fname
  # Puts quotes around each doc name and joins on commas
  doc_names_str = doc_names.map { |str| %Q{"#{str}"} }.join(", ")

  %Q{
num.topics <- function(US, VS)
{
    max.topics <- 5

    min(nrow(US), nrow(VS), max.topics)
}

biplot2 <- function(doc.scores,
                    term.loadings,
                    doc.names,
                    num.terms.to.keep = 0,
                    topic.x = 1,
                    topic.y = 2,
                    doc.color=rgb(0, 0, 0, 0.7),
                    term.color=rgb(1, 0, 0, 0.35))
{
    ## Lims for square plot that includes the origin
    all.points <- unlist(c(0,
                           doc.scores[, c(topic.x, topic.y)],
                           term.loadings[, c(topic.x, topic.y)]))
    lims <- c(min(all.points),
              max(all.points)) * 1.1

    par(mfrow=c(1,1), lwd=2, cex.axis=0.8, cex.lab=1.2)

    plot(doc.scores[, c(topic.x, topic.y)],
         xlab="",
         ylab="",
         bty="n",
         axes=F,
         type="n",
         xlim=lims,
         ylim=lims)

    grid(lwd=1)
    box()
    axis(1)
    axis(2, las=1)
    title(xlab = paste("Topic", topic.x, sep=" "),
          ylab = paste("Topic", topic.y, sep=" "),
          line=2.5)

    if (num.terms.to.keep == 0) {
        n <- nrow(term.loadings)
    } else {
        n <- round(num.terms.to.keep / 2)
    }

    abs.term.loadings <- apply(term.loadings, MARGIN=2, FUN=abs)
    top.term.loadings.x <- order(abs.term.loadings[, topic.x], decreasing=T)[1:n]
    top.term.loadings.y <- order(abs.term.loadings[, topic.y], decreasing=T)[1:n]

    top.term.loadings <- rbind(term.loadings[top.term.loadings.x, ],
                               term.loadings[top.term.loadings.y, ])

    arrows(0,
           0,
           x1=top.term.loadings[, topic.x],
           y1=top.term.loadings[, topic.y],
           col=term.color,
           length=0.08)

    offset <- abs(lims[1] - lims[2]) / 35

    actual.scores <- cbind(doc.scores[, topic.x],
                           doc.scores[, topic.y])

    offset.scores <- cbind(doc.scores[, topic.x],
                           doc.scores[, topic.y] - offset)

    points(actual.scores,
           pch=16,
           cex=0.8,
           col=doc.color)
    text(offset.scores,
         labels=doc.names,
         cex=0.8,
         col=doc.color)
}

US <- read.table("#{svd_US_fname}", header=F, sep=" ", skip=1)

VS <- read.table("#{svd_VS_fname}", header=F, sep=" ", skip=1)

n <- nrow(VS)

## Projections scaled to unit variance
doc.projections <- VS / sqrt(n - 1)

term.loadings <- US / sqrt(n - 1)

pdf("#{pdf_fname}", width=8, height=8)
num <- num.topics(US, VS)
for (x in 1:(num - 1)) {
    for (y in (x+1):num) {
        biplot2(doc.projections,
                term.loadings,
                c(#{doc_names_str}),
                num.terms.to.keep=50,
                topic.x=x,
                topic.y=y)
    }
}
invisible(dev.off())
}
end

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
  PIPELINE_VERSION = "0.9.0-alpha"
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

def run_process_svd process_svd,
                    svd_U_fname,
                    svd_S_fname,
                    svd_V_fname,
                    sing_val_inflection_point,
                    redsvd_dir,
                    log_out_fname,
                    log_err_fname

  cmd = "#{process_svd} " +
        "#{svd_U_fname} " +
        "#{svd_S_fname} " +
        "#{svd_V_fname} " +
        "#{sing_val_inflection_point} " +
        "#{redsvd_dir} " +
        "1>> #{log_out_fname} 2>> #{log_err_fname}"
  Process.run_and_time_it! "Process SVD", cmd
end


def run_td_matrix_counts td_matrix,
                         cluster_outfname,
                         outdir,
                         tf_func,
                         idf_func,
                         singleton_weight,
                         log_out_fname,
                         log_err_fname,
                         td_matrix_outfiles

  if all_files_exist? td_matrix_outfiles
    AbortIf.logger.info { "Term-doc matrix already created, skipping" }
  else
    cmd = "#{td_matrix} " +
          "#{cluster_outfname} " +
          "#{outdir} " +
          "#{tf_func} " +
          "#{idf_func} " +
          "#{singleton_weight} " +
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

# Returns the (index + 1) of the inflection point. Doubles as number
# of items until the inflection point.
def inflection_point dat
  if dat.length <= 3
    dat.length
  else
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
end

THIS_DIR = File.absolute_path(File.join File.dirname __FILE__)

opts = Trollop.options do
  version Lsa::VERSION_BANNER

  banner <<-EOS

#{Lsa::VERSION_BANNER}

  Latent semantic analysis pipeline for genomes and metagenomes.

  --num-topics is the number of topics to use for the distance
    calculations. Pass 0 for automatic mode. If --num-topics is
    greater than the maximum possible topics, the maximum number of
    topics is used instead.

  --percent-of-terms-per-topic has an automatic mode, pass 0

  tf_func options: tf_raw, tf_binary, tf_freq, tf_log_norm

  idf_func options: idf, idf_smooth, idf_const

  true singleton weights
    - pass 0 to completely drop true singletons
    - any value between 0 and 1: decrease their contribution (try 0.5
      to reduce their weight by half)
    - pass 1 to let them be (this is default)
    - any value greater than 1: increase their contribution

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
      "The maximum number of topics to use for distance calculation" +
      " (Use 0 for automatic.)",
      default: 0)
  opt(:percent_of_terms_per_topic,
      "What percentage of top terms per topic do you want to look at?" +
      " (Use 0 for automatic)",
      default: 0)

  opt(:tf_func,
      "Term frequency function",
      default: "tf_raw")
  opt(:idf_func,
      "Inverse document frequency function",
      default: "idf_smooth")
  opt(:singleton_weight,
      "Weight true singletons by this value",
      default: 1.0,
      type: :double)
end

tf_func_opts = ["tf_raw",
                "tf_binary",
                "tf_freq",
                "tf_log_norm"]
idf_func_opts = ["idf",
                 "idf_smooth",
                 "idf_const"]

abort_unless tf_func_opts.include?(opts[:tf_func]),
             "--tf-func must be one of #{tf_func_opts}. Got '#{opts[:tf_func]}'"
abort_unless idf_func_opts.include?(opts[:idf_func]),
             "--tf-func must be one of #{idf_func_opts}. Got '#{opts[:idf_func]}'"
abort_unless opts[:singleton_weight] >= 0,
             "--singleton-weight must be >= 0. Got #{opts[:singleton_weight]}"

abort_unless opts[:percent_of_terms_per_topic] >= 0 &&
             opts[:percent_of_terms_per_topic] <= 100,
             "--percent-of-terms-per-topic must be >= 0 and <= 100"

abort_unless Dir.exist?(opts[:binary_dir]),
             "The directory specified by --bin-dir doesn't exist"

abort_unless opts[:num_topics] >= 0,
             "--num-topics must be >= 0"

######################################################################
# check commands
################


prep_seq_files = File.join opts[:binary_dir], "prep_seq_files.rb"
cluster = File.join opts[:binary_dir], "cluster.rb"
td_matrix = File.join opts[:binary_dir], "td_matrix"
redsvd = File.join opts[:binary_dir], "redsvd"
mmseqs = opts[:mmseqs]
make_color_maps = File.join opts[:binary_dir], "make_color_maps.rb"
process_svd = File.join opts[:binary_dir], "process_svd"

abort_unless_command prep_seq_files
abort_unless_command cluster
abort_unless_command td_matrix
abort_unless_command redsvd
abort_unless_command mmseqs
abort_unless_command make_color_maps
abort_unless_command process_svd

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

  biplot_rscript_fname =
    File.join rscript_dir, "biplots.r"
  biplot_pdf_fname =
    File.join redsvd_dir, "biplots.pdf"

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
                       opts[:tf_func],
                       opts[:idf_func],
                       opts[:singleton_weight],
                       log_out_fname,
                       log_err_fname,
                       td_matrix_outfiles

  # If we are in the original metadata group, move the mapping files
  # into the trees dir
  if opts[:mapping] && metadata_group_label == "original"
    # TODO doesn't feel right
    if Dir.exist? trees_dir
      FileUtils.mv Dir.glob(File.join(color_map_dir, "*")), trees_dir
    else
      FileUtils.mv color_map_dir, trees_dir
    end
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
  # if opts[:num_topics].zero?
  #   opts[:num_topics] = max_topics num_terms, num_docs
  # end

  # LSA
  # Always calculate the max number of topics in the LSA step
  run_redsvd redsvd,
             td_matrix_outfname,
             redsvd_outf_base,
             max_topics(num_terms, num_docs),
             log_out_fname,
             log_err_fname

  sing_vals = File.open(svd_S_fname, "rt").
              read.
              chomp.
              split("\n").
              drop(1). # The first line is the number of sing_vals
              map(&:to_f)

  if opts[:num_topics].zero? # automatic mode
    topics_for_dist_calc = inflection_point sing_vals
  elsif opts[:num_topics] > max_topics(num_terms, num_docs)
    AbortIf.logger.warn { "--num-topics #{opts[:num_topics]} was > than max topics (#{max_topics(num_terms, num_docs)}). Using max topics" }
    topics_for_dist_calc = max_topics(num_terms, num_docs)
  else
    topics_for_dist_calc = opts[:num_topics]
  end

  run_process_svd process_svd,
                  svd_U_fname,
                  svd_S_fname,
                  svd_V_fname,
                  topics_for_dist_calc,
                  redsvd_dir,
                  log_out_fname,
                  log_err_fname

  # Make biplots
  Time.time_it "Making biplots", AbortIf.logger do
    File.open(biplot_rscript_fname, "w") do |f|
      doc_names = idx2doc.
                  sort_by { |idx, name| idx }.
                  map { |idx, name| name }

      f.puts biplot_rscript svd_US_fname,
                            svd_VS_fname,
                            doc_names,
                            biplot_pdf_fname
    end

    Process.run_and_time_it! "Plotting biplots",
                             "Rscript #{biplot_rscript_fname}"
  end

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
           drop(1). # Line specifying nrows, ncols
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

      # Remove the doc from the term name
      # TODO this will break if headers are not unique across original files?
      terms_no_doc_annotation =
        top_terms.map(&:first).map { |term| term.split("~").last }

      # Keep the term names in the ht
      topic2top_terms[topic_idx] = Set.new terms_no_doc_annotation

      top_terms.each do |term, weight, abs_weight|
        f.puts [topic_idx, term, weight].join " "
      end
    end
  end

  # Plot the term weights per topic plots
  rscript_str = %Q{
#{PLOT_FUNCTION}

dat <- read.table("#{svd_US_fname}", sep=" ", skip=1)

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
    File.open(svd_VS_US_dis_fname, "rt").each_line.with_index do |line, line_idx|
      unless line_idx.zero? # skip the header line
        doc_idx = line_idx - 1
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
        cluster_centers_fname =
          File.join doc_dir,
                    "term_doc_dists_cluster_centers.txt"

        kmeans_r_fname =
          File.join rscript_dir,
                    "term_doc_dists_kmeans.doc_#{doc_idx}_#{doc_name}.r"
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

        rscript_str = %Q{
#{PLOT_FUNCTION}

## num.centroids <- #{num_centroids}
dat <- read.table("#{dists_fname}", col.names=c("dist"))
dat.sorted <- dat[with(dat, order(dist)),]
dist.counts <- c(#{dist_counts.map { |dist, count| count }.join(", ")})

## Merge centroids that are really close
centroids <- c(#{centroids.join(", ")})
tolerance <- 2.5 / 100
hc <- hclust(dist(centroids))
groups <- cutree(hc, h = tolerance)
centroid.groups <- as.data.frame(cbind(centroids, groups))
merged.centroids <- unlist(lapply(sort(unique(centroid.groups$groups)), function(group) {
    centroid.set <- subset(centroid.groups, centroid.groups$groups == group)$centroids
    mean(centroid.set)
}))
num.centroids <- length(merged.centroids)

## Get the clusters
k <- kmeans(dat, merged.centroids)
write.table(k$cluster, "#{clusters_fname}", row.names=F, col.names=F, quote=F)
write.table(k$centers, "#{cluster_centers_fname}", row.names=F, col.names=F, quote=F)

dat.with.centroid <- as.data.frame(cbind(dat, kcluster=k$cluster))

## Plot dist counts to show how number of centroids was determined
## TODO should make this more tolerant like call a mean dist equal if
## within a certain tolerance?
pdf("#{centroids_plot_fname}", width=8, height=5)
plot.colored.by.inflection.point(dist.counts, num.centroids, ylab="Distance")
invisible(dev.off())

## Plot how the data fits into the clusters
pdf("#{dists_plot_fname}", width=8, height=5)
par(mfrow=c(1,1), lwd=2, cex.axis=0.8, cex.lab=1.2)
dat.with.centroid.sorted <- dat.with.centroid[order(dat.with.centroid$dist), ]
plot(1:nrow(dat.with.centroid.sorted), dat.with.centroid.sorted$dist, xlab="", ylab="", bty="n", axes=F, type="n")
grid(lwd=1)
box()
axis(1)
axis(2, las=1)
title(xlab="Rank", line=2.5)
title(ylab="Count", line=2.75)
points(1:nrow(dat.with.centroid.sorted), dat.with.centroid.sorted$dist, col=rainbow(num.centroids)[dat.with.centroid.sorted$kcluster], pch=16, cex=0.8)
invisible(dev.off())
}

        File.open(kmeans_r_fname, "w") do |f|
          f.puts rscript_str
        end

        Process.run_and_time_it! "Running kmeans", "Rscript #{kmeans_r_fname}"

        # This is 1-based to match the R clusters
        cluster2center = File.open(cluster_centers_fname, "rt").
                         read.
                         chomp.
                         split("\n").
                         map.with_index { |center, idx| [idx + 1, center] }.
                         to_h

        centers = cluster2center.values

        AbortIf.logger.info { "#{doc_name} will have #{centers.count} groups for top terms by doc" }

        center2new_idx = centers.sort_by { |center| center.to_f }.
                         map.
                         with_index { |center, idx| [center, idx + 1] }.
                         to_h

        clusters = File.open(clusters_fname, "rt").
                   read.
                   chomp.
                   split("\n").
                   map { |cluster| center2new_idx[cluster2center[cluster.to_i]] }

        dists_w_cluster =
          dists.map.with_index do |dist, term_idx|
          [clusters[term_idx], dist, term_idx]
        end

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
      nrows, ncols = line.chomp.split(" ").map(&:to_i)
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

        # Make the cluster file for this doc
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
      if headers.include? header_no_tilde
        # seqs can be top seq in multiple topics
        top_terms_per_topic_fnames[topic].puts rec
      end
    end
  end

  # Close the files
  top_terms_per_topic_fnames.each do |topic, f|
    f.close
  end

  term_doc_dist_cluster_fnames.each do |doc, outfiles|
    outfiles.each do |cluster, outf|
      outf.close
    end
  end
end

##################
# run the pipeline
######################################################################

AbortIf.logger.info { "Done!" }

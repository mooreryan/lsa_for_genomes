#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require "aai"
require "abort_if"
require "fileutils"
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
  PIPELINE_VERSION = "0.4.1"
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

THIS_DIR = File.absolute_path(File.join File.dirname __FILE__)

opts = Trollop.options do
  version Lsa::VERSION_BANNER

  banner <<-EOS

#{Lsa::VERSION_BANNER}

  Latent semantic analysis pipeline for genomes and metagenomes.

  --num-topics is the LIMIT on topics. You may see fewer depending on
    the data. If you pass --num-topics 0, all topics will be used.

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
      " (Use zero for maximum number of topics.)",
      default: 0)
  opt(:percent_of_terms,
      "What percentage of top terms do you want to look at?",
      default: 1)
end

KEEP_PERCENT = opts[:percent_of_terms]
abort_unless KEEP_PERCENT > 0 && KEEP_PERCENT <= 100,
             "--how-many-terms must be > 0 and <= 100"

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

abort_unless_command prep_seq_files
abort_unless_command cluster
abort_unless_command td_matrix
abort_unless_command redsvd
abort_unless_command mmseqs

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

# Cluster

if all_files_exist? mmseqs_outfiles
  AbortIf.logger.info { "Clustering already done, skipping" }
else
  cmd = "#{cluster} #{mmseqs} #{opts[:cpus]} #{prepped_seq_files} #{cluster_dir} #{prepped_seq_files_basename} 1>> #{log_out_fname} 2>> #{log_err_fname}"
  Process.run_and_time_it! "Clustering ORFs", cmd
end

label2outf = nil
doc2new_doc = nil
new_cluster_outfnames = []
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
end

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
  tmp_r_script_fname =
    File.join rscript_dir, "make_tree.r"
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
  terms_closest_to_docs =
    File.join top_terms_by_doc_dir, "top_terms_by_doc.txt"

  run_td_matrix_counts td_matrix,
                       cluster_outfname,
                       File.dirname(td_matrix_outfname),
                       log_out_fname,
                       log_err_fname,
                       td_matrix_outfiles

  # Parse the idx to name map files
  idx2doc  = parse_index_mapping_file idx_to_doc_outfname
  idx2term = parse_index_mapping_file idx_to_term_outfname

  num_docs  = idx2doc.count
  num_terms = idx2term.count

  num_terms_to_keep = (num_terms * (KEEP_PERCENT / 100.0)).round
  AbortIf.logger.info do
    "Keeping #{num_terms_to_keep} of #{num_terms} terms for each doc"
  end

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

  topic2top_terms = {}

  File.open(top_terms_by_topic, "w") do |f|
    topic2top_terms = {}
    all_weights = []
    File.open(svd_U_fname, "rt").each_line.with_index do |line, term_idx|
      *weights = line.chomp.split(" ").map { |weight| weight.to_f }

      weights_with_term_idx = weights.map { |weight| [term_idx, weight, weight.abs] }

      all_weights << weights_with_term_idx
      # TODO assert correct number of topics
    end

    num_topics = all_weights.first.count
    num_topics.times do |topic_idx|
      topic2top_terms[topic_idx] = Set.new
    end

    # Now get the top terms for the topics. In this case, the higher
    # the weight, the better it is.
    all_weights.
      transpose. # iterate over topics
      map do |elem|
      elem.
        sort_by { |term_idx, weight, abs_weight| abs_weight }.
        reverse.
        take(num_terms_to_keep) # we want only the num_terms_to_keep highest weights per topic
    end.each_with_index do |elem, topic_index| # each element holds weights for topics
      elem.each do |term_idx, weight, abs_weight|
        abort_unless idx2term.has_key?(term_idx),
                     "#{term_idx} missing from idx2term"

        term_name = idx2term[term_idx]
        topic2top_terms[topic_index] << term_name

        f.puts [topic_index, term_name, weight].join " "
      end
    end
  end


  AbortIf.logger.info { "Finding top terms for each document" }
  doc2top_terms = {}

  File.open(terms_closest_to_docs, "w") do |f|
    # Get genes closest to each doc
    File.open(svd_VS_US_dis_fname, "rt").each_line.with_index do |line, doc_idx|
      *dists = line.chomp.split " "

      abort_unless idx2doc.has_key?(doc_idx),
                   "Doc index #{doc_idx} is missing from the idx2doc hash table"

      doc_name = idx2doc[doc_idx]

      abort_unless dists.count == idx2term.count,
                   "Number of terms from #{svd_VS_US_dis_fname} does not match number of terms in idx2term hash table."

      dists_with_idx =
        dists.map.with_index { |dist, term_idx| [term_idx, dist] }

      # lowest distances are the best
      closest_terms =
        dists_with_idx.sort_by { |term_idx, dist| dist }.
        take(num_terms_to_keep)

      doc2top_terms[doc_name] = Set.new
      closest_terms.each do |term_idx, dist|
        abort_unless idx2term.has_key?(term_idx),
                     "Missing term index #{term_idx} from the idx2term hash table."

        term_name = idx2term[term_idx]
        doc2top_terms[doc_name] << term_name
        f.puts [doc_name, term_name, dist].join " "
      end
    end
  end

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

print("Reading doc dist")
proj.docs.dist <- as.dist(read.table("#{svd_VS_dis_fname_for_r}", header=T, sep=" "))

print("Making doc tree")
proj.docs.dist.tree <- as.phylo(hclust(proj.docs.dist, method="average"))

print("Printing doc tree")
write.tree(proj.docs.dist.tree, file="#{newick_docs_fname}")

}

  File.open(tmp_r_script_fname, "w") do |f|
    f.puts rscript_str
  end

  cmd = "Rscript #{tmp_r_script_fname}"
  Process.run_and_time_it! "Generating trees", cmd

  top_terms_per_topic_fnames = {}
  topic2top_terms.each do |topic, terms|
    fname = File.join top_terms_by_topic_dir,
                      "top_terms.topic_#{topic.to_i}.fa"

    top_terms_per_topic_fnames[topic] = File.open fname, "w"
  end

  top_terms_per_doc_fnames = {}
  doc2top_terms.each do |doc, terms|
    fname = File.join top_terms_by_doc_dir,
                      "top_terms.doc_#{doc}.fa"

    top_terms_per_doc_fnames[doc] = File.open fname, "w"
  end

  AbortIf.logger.info { "Grepping seq ids" }
  # Grep them from the prepped file
  ParseFasta::SeqFile.open(prepped_seq_files).each_record do |rec|
    doc2top_terms.each do |doc, headers|
      if headers.include? rec.header
        top_terms_per_doc_fnames[doc].puts rec
      end
    end

    topic2top_terms.each do |topic, headers|
      if headers.include? rec.header
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

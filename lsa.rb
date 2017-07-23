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
  PIPELINE_VERSION   = "0.3.0"
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

THIS_DIR = File.join File.dirname __FILE__

opts = Trollop.options do
  version Lsa::VERSION_BANNER

  banner <<-EOS

#{Lsa::VERSION_BANNER}

  Latent semantic analysis pipeline for genomes and metagenomes.

  --num-topics is the LIMIT on topics. You may see fewer depending on
    the data.

  Options:
  EOS

  opt(:bin_dir,
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
      default: File.join(THIS_DIR, "lsa_output"))
  opt(:mapping,
      "Mapping file (optional)",
      type: :string)

  opt(:cpus,
      "Number of CPUs to use",
      default: 3)

  opt(:num_topics,
      "The maximum number of topics to calculate for LSA",
      default: 20)
end

abort_unless Dir.exist?(opts[:bin_dir]),
             "The directory specified by --bin-dir doesn't exist"

######################################################################
# check commands
################


prep_seq_files = File.join opts[:bin_dir], "prep_seq_files.rb"
cluster = File.join opts[:bin_dir], "cluster.rb"
td_matrix = File.join opts[:bin_dir], "td_matrix"
lsa_py = File.join opts[:bin_dir], "lsa.py"
mmseqs = opts[:mmseqs]

abort_unless_command prep_seq_files
abort_unless_command cluster
abort_unless_command td_matrix
abort_unless_command lsa_py
abort_unless_command mmseqs

################
# check commands
######################################################################

######################################################################
# outfiles
##########

prepped_seq_files = File.join opts[:outdir], "all_prepped.fa"

mmseqs_final_outf = "#{prepped_seq_files}.clu.tsv.sorted"
mmseqs_outfiles = [
  "#{prepped_seq_files}.DB",
  "#{prepped_seq_files}.DB.index",
  "#{prepped_seq_files}.DB.lookup",
  "#{prepped_seq_files}.DB_h",
  "#{prepped_seq_files}.DB_h.index",
  "#{prepped_seq_files}.clu",
  "#{prepped_seq_files}.clu.index",
  "#{prepped_seq_files}.clu.tsv",
  mmseqs_final_outf,
]


outf = File.join opts[:outdir], "lsa_out.txt"
errf = File.join opts[:outdir], "lsa_err.txt"

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
  cmd = "#{prep_seq_files} #{infiles} 1> #{prepped_seq_files} 2> #{errf}"
  Process.run_and_time_it! "Prepping ORFs", cmd
  abort_unless_file_exists prepped_seq_files
end

# Cluster

if all_files_exist? mmseqs_outfiles
  AbortIf.logger.info { "Clustering already done, skipping" }
else
  cmd = "#{cluster} #{mmseqs} #{opts[:cpus]} #{prepped_seq_files} 1>> #{outf} 2>> #{errf}"
  Process.run_and_time_it! "Clustering ORFs", cmd
end

label2outf = nil
doc2new_doc = nil
new_cluster_outfnames = []
if opts[:mapping]
  # Read the mapping file
  AbortIf.logger.info { "Parsing mapping file" }
  label2outf, doc2new_doc = Lsa::parse_mapping_file opts[:mapping],
                                                    mmseqs_final_outf

  # Make the new cluster files
  AbortIf.logger.info { "Making new cluster files" }
  new_cluster_outfnames = Lsa::make_new_cluster_files mmseqs_final_outf,
                                                      label2outf,
                                                      doc2new_doc
end

# From here on, run everything in a big loop for all the cluster files
[mmseqs_final_outf, new_cluster_outfnames].flatten.each do |cluster_outf|
  AbortIf.logger.info { "Running big loop on #{cluster_outf}" }


  td_matrix_outf = "#{cluster_outf}.seanie_lsa.td_matrix.mm.txt"
  idx_to_doc_outf = "#{cluster_outf}.seanie_lsa.idx_to_doc_map.txt"
  idx_to_term_outf = "#{cluster_outf}.seanie_lsa.idx_to_term_map.txt"
  td_matrix_outfiles = [
    td_matrix_outf,
    idx_to_doc_outf,
    idx_to_term_outf,
  ]

  lsa_py_dist_outf = "#{cluster_outf}.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt"
  lsa_py_outfiles = [
    "#{cluster_outf}.seanie_lsa.td_matrix.mm.txt.lsa_py.rows_are_terms.txt",
    "#{cluster_outf}.seanie_lsa.td_matrix.mm.txt.lsa_py.singular_values.txt",
    "#{cluster_outf}.seanie_lsa.td_matrix.mm.txt.lsa_py.top_terms.txt",
    lsa_py_dist_outf,
    "#{cluster_outf}.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_matrix.txt",
  ]

  tmp_r_script_fname = File.join "#{cluster_outf}.make_tree.r"
  newick_outf = lsa_py_dist_outf + ".newick.txt"

  # Tf-idf counts

  if all_files_exist? td_matrix_outfiles
    AbortIf.logger.info { "Term-doc matrix already created, skipping" }
  else
    cmd = "#{td_matrix} #{cluster_outf} 1>> #{outf} 2>> #{errf}"
    Process.run_and_time_it! "Building term-doc matrix", cmd
  end

  # LSA

  if all_files_exist? lsa_py_outfiles
    AbortIf.logger.info { "LSA tranform already done, skipping" }
  else
    cmd = "#{lsa_py} #{opts[:num_topics]} #{td_matrix_outf} #{idx_to_term_outf} 1>> #{outf} 2>> #{errf}"
    Process.run_and_time_it! "Running the LSA transform", cmd
  end

  # Make the R script for the tree

  idx2doc = {}
  File.open(idx_to_doc_outf, "rt").each_line do |line|
    idx, doc = line.chomp.split "\t"

    idx2doc[idx.to_i] = doc
  end

  col_names = idx2doc.
              sort_by { |idx, doc| idx }.
              map { |idx, doc| "\"#{doc}\"" }.
              join ", "

  rscript_str = %Q{
library("ape")

transformed.doc.dist <- as.dist(read.table("#{lsa_py_dist_outf}", col.names=c(#{col_names})))
transformed.doc.dist.tree <- as.phylo(hclust(transformed.doc.dist, method="average"))

write.tree(transformed.doc.dist.tree, file="#{newick_outf}")
}

  File.open(tmp_r_script_fname, "w") do |f|
    f.puts rscript_str
  end

  cmd = "Rscript #{tmp_r_script_fname}"
  Process.run_and_time_it! "Making tree", cmd

  # Get the actual sequences for top seqs
  top_terms_fname = "#{td_matrix_outf}.lsa_py.top_terms.txt"
  ids_to_grep_fname = "#{top_terms_fname}.ids.txt"

  AbortIf.logger.info { "Reading top terms" }
  topic2seqid = {}
  File.open(top_terms_fname, "rt").each_line do |line|
    unless line.start_with? "topic"
      topic, term, val = line.chomp.split " "

      if topic2seqid.has_key? topic
        topic2seqid[topic] << term
      else
        topic2seqid[topic] = Set.new [term]
      end
    end
  end

  # Make the outfiles
  top_fasta_fnames = {}
  topic2seqid.each do |topic, seqid|
    top_fasta_fnames[topic] =
      File.open("#{ids_to_grep_fname}.topic_#{topic.to_i+1}_seqs.fa", "w")
  end

  AbortIf.logger.info { "Grepping seq ids" }
  # Grep them from the prepped file
  ParseFasta::SeqFile.open(prepped_seq_files).each_record do |rec|
    topic2seqid.each do |topic, headers|
      if headers.include? rec.header
        # seqs can be top in multiple topics
        top_fasta_fnames[topic].puts rec
      end
    end
  end

  # Close the files
  top_fasta_fnames.each do |topic, f|
    f.close
  end
end

##################
# run the pipeline
######################################################################

AbortIf.logger.info { "Done!" }

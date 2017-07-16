#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require "aai"
require "abort_if"
require "fileutils"
require "trollop"

include AbortIf
Process.extend Aai::CoreExtensions::Process

module Aai
  extend Aai
  extend Aai::Utils
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

THIS_DIR = File.join File.dirname __FILE__

opts = Trollop.options do
  banner <<-EOS

  Latent semantic analysis pipeline for genomes and metagenomes

  Options:
  EOS

  opt(:bin_dir,
      "Folder with the LSA scripts and binaries",
      default: File.join(THIS_DIR, "bin"))

  opt(:infiles,
      "Files with ORF clusters",
      type: :strings)
  opt(:outdir,
      "Output directory",
      type: :string,
      default: File.join(THIS_DIR, "lsa_output"))

  opt(:cpus,
      "Number of CPUs to use",
      default: 3)

  # opt(:force,
  #     "Overwrite contents of outdir",
  #     default: false)
end

######################################################################
# check commands
################


prep_seq_files = File.join opts[:bin_dir], "prep_seq_files.rb"
cluster = File.join opts[:bin_dir], "cluster.rb"
td_matrix = File.join opts[:bin_dir], "td_matrix"
lsa_py = File.join opts[:bin_dir], "lsa.py"

abort_unless_command prep_seq_files
abort_unless_command cluster
abort_unless_command td_matrix
abort_unless_command lsa_py

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

td_matrix_outf = "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt"
idx_to_doc_outf = "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.idx_to_doc_map.txt"
idx_to_term_outf = "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.idx_to_term_map.txt"
td_matrix_outfiles = [
  td_matrix_outf,
  idx_to_doc_outf,
  idx_to_term_outf,
]

lsa_py_outfiles = [
  "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.rows_are_terms.txt",
  "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.singular_values.txt",
  "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.top_terms.txt",
  "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt",
  "#{prepped_seq_files}.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_matrix.txt",
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
  cmd = "#{cluster} #{opts[:cpus]} #{prepped_seq_files} 1>> #{outf} 2>> #{errf}"
  Process.run_and_time_it! "Clustering ORFs", cmd
end

# Tf-idf counts

if all_files_exist? td_matrix_outfiles
  AbortIf.logger.info { "Term-doc matrix already created, skipping" }
else
  cmd = "#{td_matrix} #{mmseqs_final_outf} 1>> #{outf} 2>> #{errf}"
  Process.run_and_time_it! "Building term-doc matrix", cmd
end

# LSA

if all_files_exist? lsa_py_outfiles
  AbortIf.logger.info { "LSA tranform already done, skipping" }
else
  cmd = "#{lsa_py} #{td_matrix_outf} #{idx_to_term_outf} 1>> #{outf} 2>> #{errf}"
  Process.run_and_time_it! "Running the LSA transform", cmd
end

##################
# run the pipeline
######################################################################

AbortIf.logger.info { "Done!" }

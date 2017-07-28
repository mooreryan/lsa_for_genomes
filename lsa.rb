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
      default: 4)

  opt(:num_topics,
      "The maximum number of topics to calculate for LSA",
      default: 20)
  opt(:percent_of_terms,
      "What percentage of top terms do you want to look at?",
      default: 1)
end

KEEP_PERCENT = opts[:percent_of_terms]
abort_unless KEEP_PERCENT > 0 && KEEP_PERCENT <= 100,
             "--how-many-terms must be > 0 and <= 100"

abort_unless Dir.exist?(opts[:bin_dir]),
             "The directory specified by --bin-dir doesn't exist"

######################################################################
# check commands
################


prep_seq_files = File.join opts[:bin_dir], "prep_seq_files.rb"
cluster = File.join opts[:bin_dir], "cluster.rb"
td_matrix = File.join opts[:bin_dir], "td_matrix"
redsvd = File.join opts[:bin_dir], "redsvd"
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


  td_matrix_outf = "#{cluster_outf}.seanie_lsa.td_matrix.txt"
  idx_to_doc_outf = "#{cluster_outf}.seanie_lsa.idx_to_doc_map.txt"
  idx_to_term_outf = "#{cluster_outf}.seanie_lsa.idx_to_term_map.txt"
  td_matrix_outfiles = [
    td_matrix_outf,
    idx_to_doc_outf,
    idx_to_term_outf,
  ]

  redsvd_outf_base = "#{td_matrix_outf}.svd"
  tmp_r_script_fname = File.join "#{cluster_outf}.make_tree.r"
  newick_terms_fname = redsvd_outf_base + ".terms_newick.txt"
  newick_docs_fname = redsvd_outf_base + ".docs_newick.txt"

  # Tf-idf counts

  if all_files_exist? td_matrix_outfiles
    AbortIf.logger.info { "Term-doc matrix already created, skipping" }
  else
    cmd = "#{td_matrix} #{cluster_outf} 1>> #{outf} 2>> #{errf}"
    Process.run_and_time_it! "Building term-doc matrix", cmd
  end



  # Parse the idx to name map files
  idx2doc = {}
  File.open(idx_to_doc_outf, "rt").each_line do |line|
    idx, doc = line.chomp.split "\t"

    idx2doc[idx.to_i] = doc
  end

  doc_names = idx2doc.
              sort_by { |idx, doc| idx }.
              map { |idx, doc| "\"#{doc}\"" }.
              join "\n"

  ordered_doc_names = idx_to_doc_outf + ".for_r"
  File.open(ordered_doc_names, "w") do |f|
    f.puts doc_names
  end

  idx2term = {}
  File.open(idx_to_term_outf, "rt").each_line do |line|
    idx, term = line.chomp.split "\t"

    idx2term[idx.to_i] = term
  end

  term_names = idx2term.
              sort_by { |idx, term| idx }.
              map { |idx, term| "\"#{term}\"" }.
              join "\n"

  ordered_term_names = idx_to_term_outf + ".for_r"
  File.open(ordered_term_names, "w") do |f|
    f.puts term_names
  end




  # LSA

  # if all_files_exist? lsa_py_outfiles
  #   AbortIf.logger.info { "LSA tranform already done, skipping" }
  # else
  cmd = "#{redsvd} --format sparse --method SVD --input #{td_matrix_outf} --output #{redsvd_outf_base} --rank #{opts[:num_topics]} 1>> #{outf} 2>> #{errf}"
  Process.run_and_time_it! "Running the LSA transform", cmd
  # end

  # Make the R script for the tree

  svd_U_fname = "#{redsvd_outf_base}.U"
  svd_S_fname = "#{redsvd_outf_base}.S"
  svd_V_fname = "#{redsvd_outf_base}.V"
  svd_VS_fname = "#{redsvd_outf_base}.VS"
  svd_US_fname = "#{redsvd_outf_base}.US"
  svd_VS_dis_fname = "#{redsvd_outf_base}.VS.dis"
  svd_VS_dis_fname_for_r = "#{redsvd_outf_base}.VS.dis.for_r"
  svd_VS_US_dis_fname = "#{redsvd_outf_base}.VS_to_US.dis"

  num_to_keep = (idx2term.count * (KEEP_PERCENT / 100.0)).round
  AbortIf.logger.info { "Keeping #{num_to_keep} of #{idx2term.count} terms for each doc" }

  AbortIf.logger.info { "Finding top terms for each document" }
  terms_closest_to_docs = svd_VS_US_dis_fname + ".terms_closest_to_docs.txt"
  File.open(terms_closest_to_docs, "w") do |f|
    # Get genes closest to each doc
    File.open(svd_VS_US_dis_fname, "rt").each_line.with_index do |line, doc_idx|
      *dists = line.chomp.split " "

      abort_unless idx2doc.has_key?(doc_idx),
                   "Doc index #{doc_idx} is missing from the idx2doc hash table"

      doc_name = idx2doc[doc_idx]

      abort_unless dists.count == idx2term.count,
                   "Number of terms from #{svd_VS_US_dis_fname} does not match number of terms in idx2term hash table."

      dists_with_idx = dists.map.with_index { |dist, term_idx| [term_idx, dist] }

      # lowest distances are the best
      closest_terms = dists_with_idx.sort_by { |term_idx, dist| dist }.take num_to_keep

      closest_terms.each do |term_idx, dist|
        abort_unless idx2term.has_key?(term_idx),
                     "Missing term index #{term_idx} from the idx2term hash table."

        term_name = idx2term[term_idx]
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

  top_terms_fname = "#{td_matrix_outf}.top_terms_by_topic.txt"

  if File.exist? top_terms_fname
    AbortIf.logger.info { "#{top_terms_fname} already exists. Deleting it." }

    File.delete top_terms_fname
  end

  rscript_str = %Q{
library("lsa")
library("ape")

write.top.terms <- function(dat, file, n=10)
{
    if (n > nrow(dat)) {
        stop("n must be <= nrow(dat)")
    }

    abs.dat <- abs(dat)

    rows <- c()
    for (topic in 1:ncol(dat)) {
        print(topic)
        top.n <- (abs.dat[order(abs.dat[, topic],
                                decreasing=T), ][1:n, ])
        seq.names <- paste(rownames(top.n), collapse="\\t")

        write(paste(topic, seq.names, sep="\\t"),
              file=file,
              ## sep="\t",
              append=T)
    }
}

print("Reading idx2term")
term.names <- as.vector(read.table("#{ordered_term_names}", sep=" ", col.names=c("terms"))$terms)

print("Reading U")
## Rows are terms, columns are the eigenvectors
U <- as.matrix(read.table("#{svd_U_fname}", sep=" ", row.names=term.names))

print("Writing top terms")
## These are top terms with respect to each topic
write.top.terms(U, "#{top_terms_fname}", n=#{opts[:how_many_terms]})

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
  Process.run_and_time_it! "Generating trees and top terms", cmd

  # Get the actual sequences for top seqs
  ids_to_grep_fname = "#{top_terms_fname}.ids.txt"

  AbortIf.logger.info { "Reading top terms" }
  topic2seqid = {}
  File.open(top_terms_fname, "rt").each_line do |line|
    unless line.start_with? "topic"
      topic, *terms = line.chomp.split " "

      abort_if topic2seqid.has_key?(topic),
               "Topic #{topic} is repeated in #{top_terms_fname}"

      topic2seqid[topic] = Set.new terms
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
        # seqs can be top seq in multiple topics
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

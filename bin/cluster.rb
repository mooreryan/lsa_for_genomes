#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

require "aai"
require "abort_if"
require "fileutils"

include AbortIf
include Aai::CoreExtensions::Process

abort_unless ARGV.count == 5,
             "USAGE: ruby #{__FILE__} /path/to/mmseqs " +
             "num_threads seqs.fa outdir outbase"

# Returns a tmp dir that does not already exist
def make_tmp_dir dir
  tmp_dir = File.join dir, "lsa_cluster_tmp.#{rand}"

  # Make sure tmp_dir doesn't already exist
  while Dir.exist? tmp_dir
    tmp_dir = File.join dir, "lsa_cluster_tmp.#{rand}"
  end

  tmp_dir
end

mmseqs  = ARGV[0]
threads = ARGV[1]
fasta   = ARGV[2]
outdir  = ARGV[3]
outbase = ARGV[4]

FileUtils.mkdir_p outdir
tmp_dir = make_tmp_dir

outbase_with_dir = File.join outdir, outbase

cmd = "#{mmseqs} createdb #{fasta} #{outbase_with_dir}.DB"
Process.run_and_time_it! "Creating DB", cmd

cmd = "#{mmseqs} cluster " +
      "#{outbase_with_dir}.DB #{outbase_with_dir}.clu #{tmp_dir} " +
      "--threads #{threads} --cascaded --min-seq-id 0.3 -c 0.7"
Process.run_and_time_it! "Clustering", cmd

cmd = "#{mmseqs} createtsv " +
      "#{outbase_with_dir}.DB #{outbase_with_dir}.DB " +
      "#{outbase_with_dir}.clu #{outbase_with_dir}.clu.tsv"
Process.run_and_time_it! "Creating tsv", cmd

# -k1,1 is to fix the bug where the svd file sometimes has more terms
# -than actual terms
cmd = "sort -k1,1 #{outbase_with_dir}.clu.tsv " +
      "> #{outbase_with_dir}.clu.tsv.sorted"
Process.run_and_time_it! "Sort tsv file", cmd

Process.run_it! "rm -r #{tmp_dir}"

AbortIf.logger.info { "#{__FILE__} done!" }

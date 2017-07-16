#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

require "aai"
require "abort_if"

include AbortIf
include Aai::CoreExtensions::Process

abort_unless ARGV.count == 2,
             "USAGE: ruby #{__FILE__} num_threads seqs.fa"

threads = ARGV[0]
fasta = ARGV[1]
tmp_dir = "arstoienarstoienarstoienarstoienarstoeinarsotien"

cmd = "mmseqs createdb #{fasta} #{fasta}.DB"
Process.run_and_time_it! "Creating DB", cmd

Process.run_it "rm -r #{tmp_dir}"
Process.run_it! "mkdir #{tmp_dir}"

cmd = "mmseqs cluster #{fasta}.DB #{fasta}.clu #{tmp_dir} " +
      "--threads #{threads} --cascaded --min-seq-id 0.3 -c 0.7"
Process.run_and_time_it! "Clustering", cmd

cmd = "mmseqs createtsv #{fasta}.DB #{fasta}.DB #{fasta}.clu " +
      "#{fasta}.clu.tsv"
Process.run_and_time_it! "Creating tsv", cmd

cmd = "sort #{fasta}.clu.tsv > #{fasta}.clu.tsv.sorted"
Process.run_and_time_it! "Sort tsv file", cmd

Process.run_it! "rm -r #{tmp_dir}"

AbortIf.logger.info { "#{__FILE__} done!" }

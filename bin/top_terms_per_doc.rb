#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")
require "abort_if"
require "set"
require "trollop"

include AbortIf

opts = Trollop.options do
  banner <<-EOS

  Process the top_terms_per_doc.txt file output by LSA.

  Writes a file with shared terms amoung the first --num-clusters
  clusters if there are any.

  Write files with unique terms for the first --num-clusters clusters
  for each document.

  Options:
  EOS

  opt(:num_clusters,
      "Number of top clusters to look at",
      default: 1)
  opt(:infile,
      "The top_terms_by_doc_fname.txt file",
      type: :string)
end



abort_unless opts[:num_clusters] >= 1,
             "--num-clusters must be at least 1"
abort_if opts[:infile].nil?,
         "--infile is a required arg. Try #{__FILE__} --help for help."
abort_unless_file_exists opts[:infile]

top_terms_by_doc_fname = opts[:infile]
num_clusters = opts[:num_clusters]
intersection_outfname = top_terms_by_doc_fname + ".intersection.txt"

docs = {}
STDERR.puts "Reading top terms file"
File.open(top_terms_by_doc_fname, "rt").each_line.with_index do |line, idx|
  STDERR.printf("Reading line: #{idx}\r") if (idx % 10000).zero?

  doc_name, doc_idx, term_name, term_idx, cluster, dist = line.chomp.split " "

  unless idx.zero?

    if cluster.to_i <= num_clusters
      if docs.has_key? doc_name
        docs[doc_name] << term_name
      else
        docs[doc_name] = Set.new([term_name])
      end
    end
  end
end

key_cluster_shared_terms = docs.values.reduce(&:&)

STDERR.puts "There were #{key_cluster_shared_terms.count} shared terms."
if key_cluster_shared_terms.count > 0
  STDERR.puts "Writing intersection"
  File.open(intersection_outfname, "w") do |f|
    key_cluster_shared_terms.sort.each do |term|
      f.puts term
    end
  end
end

total_docs = docs.count

STDERR.puts "Writing unique terms"
docs.each_with_index do |(doc, terms), idx|
  STDERR.printf "Writing doc ##{idx+1} of #{total_docs}\r"
  outfname = top_terms_by_doc_fname + ".doc_#{doc}.unique_terms.txt"
  other_docs = Set.new(docs.keys) - [doc]

  other_terms = other_docs.map do |other_doc|
    docs[other_doc]
  end.reduce(&:|)

  unique_terms = terms - other_terms

  if unique_terms.count > 0
    File.open(outfname, "w") do |f|
      unique_terms.each do |term|
        f.puts [doc, term].join(" ")
      end
    end
  end
end

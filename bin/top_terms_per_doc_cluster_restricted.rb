#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")
require "abort_if"
require "set"
require "trollop"

include AbortIf

opts = Trollop.options do
  banner <<-EOS

  Process the top_terms_per_doc.txt file output by LSA.

  For each document cluster:

  Writes a file with shared terms amoung the first --num-clusters
  clusters if there are any.

  Write files with unique terms for the first --num-clusters clusters
  for each document.

  Options:
  EOS

  opt(:num_term_clusters,
      "Number of top clusters to look at",
      default: 1)
  opt(:top_terms_by_doc,
      "The top_terms_by_doc_fname.txt file",
      type: :string)
  opt(:document_clusters,
      "Document cluters from terms_per_node.rb script",
      type: :string)
end



abort_unless opts[:num_term_clusters] >= 1,
             "--num-clusters must be at least 1"

abort_if opts[:top_terms_by_doc].nil?,
         "--top-terms-by-doc is a required arg. Try #{__FILE__} --help for help."
abort_unless_file_exists opts[:top_terms_by_doc]

abort_if opts[:document_clusters].nil?,
         "--document-clusters is a required arg. Try #{__FILE__} --help for help."
abort_unless_file_exists opts[:document_clusters]

top_terms_by_doc_fname = opts[:top_terms_by_doc]
doc_clusters_fname = opts[:document_clusters]
num_term_clusters = opts[:num_term_clusters]

STDERR.puts "Reading document clusters"
clusters = {}
File.open(doc_clusters_fname, "rt").each_line do |line|
  cluster, size, *cluster_docs = line.chomp.split "\t"

  # TODO these can have spaces, but the ones in the top terms file don't
  clusters[cluster] = cluster_docs.map { |s| s.tr(" ", "_").gsub(/_+/, "_") }
end

# TODO this name is confusing, make it doc2terms
docs = {}
STDERR.puts "Reading top terms file"
File.open(top_terms_by_doc_fname, "rt").each_line.with_index do |line, idx|
  STDERR.printf("Reading line: #{idx}\r") if (idx % 10000).zero?

  doc_name, doc_idx, term_name, term_idx, cluster, dist = line.chomp.split " "

  unless idx.zero?

    if cluster.to_i <= num_term_clusters
      if docs.has_key? doc_name
        docs[doc_name] << term_name
      else
        docs[doc_name] = Set.new([term_name])
      end
    end
  end
end

clusters.each do |cluster, cluster_docs|
  key_cluster_shared_terms = cluster_docs.map { |doc| docs[doc] }.reduce(&:&)

  STDERR.puts "There were #{key_cluster_shared_terms.count} shared terms for cluster #{cluster}."
  if key_cluster_shared_terms.count > 0
    intersection_outfname = top_terms_by_doc_fname + ".intersection.cluster_#{cluster}.txt"

    STDERR.puts "Writing intersection"
    File.open(intersection_outfname, "w") do |f|
      key_cluster_shared_terms.sort.each do |term|
        f.puts term
      end
    end
  end
end

STDERR.puts "Writing unique terms"

clusters.each do |cluster, cluster_docs|
  STDERR.puts "Working on #{cluster}"
  total_docs = cluster_docs.count

  cluster_docs.each_with_index do |doc, idx|
    terms = docs[doc]

    STDERR.printf "Writing doc ##{idx+1} of #{total_docs}\r"
    outfname = top_terms_by_doc_fname + ".unique_terms.cluster_#{cluster}.doc_#{doc}.txt"
    other_docs = Set.new(docs.keys) - [doc]

    other_terms = other_docs.map do |other_doc|
      docs[other_doc]
    end.reduce(&:|)

    unique_terms = terms - other_terms

    # if unique_terms.count > 0
    File.open(outfname, "w") do |f|
      unique_terms.each do |term|
        f.puts [doc, term].join(" ")
      end
    end
    # end
    STDERR.puts
  end
end

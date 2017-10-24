#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")
require "abort_if"
require "set"
require "trollop"

include AbortIf

opts = Trollop.options do
  banner <<-EOS

  Process the top_terms_per_topic.txt file output by LSA.

  Writes a file with shared terms that are top terms in each topic.

  Write files with unique terms among the top terms per topic.

  Options:
  EOS

  opt(:infile,
      "The top_terms_by_topic_fname.txt file",
      type: :string)
end



abort_if opts[:infile].nil?,
         "--infile is a required arg. Try #{__FILE__} --help for help."
abort_unless_file_exists opts[:infile]

top_terms_by_topic_fname = opts[:infile]
intersection_outfname = top_terms_by_topic_fname + ".intersection.txt"

topics = {}
STDERR.puts "Reading top terms file"
File.open(top_terms_by_topic_fname, "rt").each_line.with_index do |line, idx|
  STDERR.printf("Reading line: #{idx}\r") if (idx % 10000).zero?

  # Note that terms in this file are all top terms (i.e., above the
  # inflection point.)
  topic_idx, term_idx, term_name, weight = line.chomp.split " "
  weight = weight.to_f

  if topics.has_key? topic_idx
    topics[topic_idx] << term_name
  else
    topics[topic_idx] = Set.new([term_name])
  end
end

shared_terms = topics.values.reduce(&:&)

STDERR.puts "There were #{shared_terms.count} shared terms."
if shared_terms.count > 0
  STDERR.puts "Writing intersection"
  File.open(intersection_outfname, "w") do |f|
    shared_terms.sort.each do |term|
      f.puts term
    end
  end
end

total_topics = topics.count

STDERR.puts "Writing unique terms"
topics.each_with_index do |(topic, terms), idx|
  STDERR.printf "Writing topic ##{idx+1} of #{total_topics}\r"
  outfname = top_terms_by_topic_fname + ".topic_#{topic}.unique_terms.txt"
  other_topics = Set.new(topics.keys) - [topic]

  other_terms = other_topics.map do |other_topic|
    topics[other_topic]
  end.reduce(&:|)

  unique_terms = terms - other_terms

  if unique_terms.count > 0
    File.open(outfname, "w") do |f|
      unique_terms.each do |term|
        f.puts [topic, term].join(" ")
      end
    end
  end
end

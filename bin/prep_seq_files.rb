#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require "abort_if"
require "parse_fasta"

include AbortIf

def clean_str str
  str.strip.gsub(/[^\p{Alnum}_]+/, "_").gsub(/_+/, "_")
end

ARGV.each do |fname|
  AbortIf.logger.info { "Working on #{fname}" }

  ParseFasta::SeqFile.open(fname).each_record do |rec|
    puts ">#{clean_str fname}~#{clean_str rec.header}"
    puts rec.seq
  end
end

AbortIf.logger.info { "#{__FILE__} done!" }

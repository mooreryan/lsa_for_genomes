#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

current = ""
File.open(ARGV.first).each_line.with_index do |line, idx|
  STDERR.printf "PARSING -- #{idx}\r" if (idx % 100000).zero?

  if idx > 1
    row, col, val = line.chomp.split

    if current == row
      print " #{col.to_i - 1}:#{val}"
    else
      current = row
      puts unless idx == 2
      print "#{col.to_i - 1}:#{val}"
    end
  end
end
puts

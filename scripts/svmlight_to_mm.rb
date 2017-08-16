#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

current = ""

num_lines = File.open(ARGV.first).each_line.count

num_rows = 0
num_cols = 0
total_entries = 0
lines = []

File.open(ARGV.first).each_line.with_index do |line, ridx|
  num_rows += 1
  # STDERR.printf "PARSING -- #{ridx}\r" if (ridx % 100).zero?

  line.chomp.split(" ").each do |col|
    total_entries += 1
    cidx, val = col.split(":")
    if (cidx.to_i + 1) > num_cols
      num_cols = cidx.to_i + 1
    end

    lines << [ridx.to_i + 1, cidx.to_i + 1, val].join(" ")
  end
end

puts '%%MatrixMarket matrix coordinate real general'
puts [num_rows, num_cols, total_entries].join " "
puts lines

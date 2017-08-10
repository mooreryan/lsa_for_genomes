#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")
require "abort_if"
require "fileutils"

include AbortIf

IROKI_KELLY_COLOR_COUNT = 19

mapping_fname = ARGV[0]
outdir = ARGV[1]
FileUtils.mkdir_p outdir

group_names = nil
rows = []
fnames = []
File.open(ARGV.first, "rt").each_line.with_index do |line, idx|
  ary = line.chomp.split "\t"
  if idx.zero?
    _, *group_names = ary
  else
    fname, *metadata = ary
    fnames << fname
    rows << metadata
  end
end

metadata_groups = rows.transpose

group_names.each_with_index do |group, idx|
  fname = File.join outdir, "#{group.tr(" ", "_")}.color_map.txt"
  data = metadata_groups[idx]

  if data.uniq.count <= IROKI_KELLY_COLOR_COUNT
    # we can do the kelly color scheme
    group2idx = data.uniq.map.with_index { |item, idx| [item, idx+1] }.to_h

    File.open(fname, "w") do |f|
      f.puts ["doc", group].join "\t"

      fnames.each_with_index do |fname, idx|
        new_doc = data[idx]

        abort_unless group2idx.has_key?(new_doc),
                     "#{new_doc} is missing from group2idx"

        f.puts [fname, group2idx[new_doc]].join "\t"
      end
    end
  end
end

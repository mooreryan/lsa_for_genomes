#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")
require "aai"
require "abort_if"
require "fileutils"

include AbortIf
Process.extend Aai::CoreExtensions::Process

module Aai
  extend Aai
  extend Aai::Utils
end


IROKI_KELLY_COLOR_COUNT = 19

def clean_str str
  str.strip.gsub(/[^\p{Alnum}_]+/, "_").gsub(/_+/, "_")
end

mapping_fname = ARGV[0]
outdir = ARGV[1]
FileUtils.mkdir_p outdir



group_names = nil
rows = []
fnames = []
File.open(mapping_fname, "rt").each_line.with_index do |line, idx|
  ary = line.chomp.split "\t"
  if idx.zero?
    _, *group_names = ary
  else
    fname, *metadata = ary
    # Storing the original fnames
    fnames << fname
    rows << metadata
  end
end

clean2orig = {}
fnames.each do |fname|
  clean2orig[clean_str(fname)] = fname
end

# Write a name map file where the names are processed through R's
# make.names function as this is what the tree will have.
names_fname = File.join outdir, "names.tmp.txt"
clean_to_rclean_name_map_fname = File.join outdir, "clean_to_rclean_name_map.txt"
File.open(names_fname, "w") { |f| f.puts fnames.map{|s| clean_str s } }
rscript_str = %Q{
dat <- read.table("#{names_fname}", header=F)
dat.map <- cbind(dat, make.names(dat$V1, unique=T))
colnames(dat.map) <- c("clean", "r_clean")
write.table(dat.map, "#{clean_to_rclean_name_map_fname}", quote=F, sep="\t", row.names=F)
}

rscript_fname = clean_to_rclean_name_map_fname + ".r"
File.open(rscript_fname, "w") { |f| f.puts rscript_str }
cmd = "Rscript #{rscript_fname}"

Process.run_and_time_it! "Generating R names", cmd
FileUtils.rm names_fname

# Get the R to original names
orig2r = {}
File.open(clean_to_rclean_name_map_fname, "rt").each_line.with_index do |line, idx|
  unless idx.zero?
    clean_name, r_name = line.chomp.split "\t"

    orig2r[clean2orig[clean_name]] = r_name
  end
end

# Write the name map
name_map_fname = File.join outdir, "name_map.txt"
File.open(name_map_fname, "w") do |f|
  orig2r.each do |orig_name, rname|
    f.puts [rname, orig_name].join "\t"
  end
end

metadata_groups = rows.transpose

group_names.each_with_index do |group, idx|
  out_fname = File.join outdir, "#{group.tr(" ", "_")}.color_map.txt"
  data = metadata_groups[idx]

  if data.uniq.count <= IROKI_KELLY_COLOR_COUNT
    # we can do the kelly color scheme
    group2idx = data.uniq.map.with_index { |item, idx| [item, idx+1] }.to_h

    File.open(out_fname, "w") do |f|
      f.puts ["#doc", group].join "\t"

      fnames.each_with_index do |fname, idx|
        new_doc = data[idx]

        abort_unless group2idx.has_key?(new_doc),
                     "#{new_doc} is missing from group2idx"

        f.puts [orig2r[fname], group2idx[new_doc]].join "\t"
      end
    end
  end
end

FileUtils.rm clean_to_rclean_name_map_fname
FileUtils.rm rscript_fname

#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

require "aai"
require "matrix"
require "abort_if"

include AbortIf
include AbortIf::Assert

Process.extend Aai::CoreExtensions::Process

def vector_projection a, b
  (a.dot(b) / b.dot(b).to_f) * b
end

def vector_rejection a, b
  a - vector_projection(a, b)
end

def index_of_max ary
  ary.each_with_index.max.last
end

# Returns the index of the inflection point
def inflection_point dat
  num_points = dat.length
  first_point = Vector[0, dat.first]
  last_point = Vector[num_points-1, dat.last]

  b = last_point - first_point

  dists = (1..num_points-2).map do |idx|
    this_point = Vector[idx, dat[idx]]
    a = this_point - first_point
    vector_rejection(a, b).norm
  end

  index_of_max(dists) + 1 # add one for the 1.. range
end

abort_unless ARGV.count == 1,
             "USAGE: #{__FILE__} text_file_with_one_column_of_numbers"

singular_values = File.open(ARGV.first, "rt").read.chomp.split("\n").map { |str| str.chomp.to_f }.sort.reverse

ipoint = inflection_point singular_values

AbortIf.logger.info { "Number of singular values: #{singular_values.count}" }
AbortIf.logger.info { "Inflection point: #{ipoint}" }

var_plot_pdf = File.join File.dirname(ARGV.first), "sv.var_plot.pdf"

rscript_str = %Q{
## This cutoff index is 1-based
plot.colored.by.inflection.point <- function(dat, cutoff, xlab="Rank", ylab="Weight")
{
    par(mfrow=c(1,1), lwd=2, cex.axis=0.8, cex.lab=1.2)

    col1 <- rgb(1, 0, 0, 0.5)
    col2 <- rgb(0, 0, 0, 0.5)

    plot(dat,
         xlab="", ylab="",
         bty="n",
         axes=F,
         type="n")

    grid(lwd=1)
    box()
    axis(1)
    axis(2, las=1)
    title(xlab=xlab, line=2.5)
    title(ylab=ylab, line=2.75)

    points(x=1:cutoff,
           y=dat[1:cutoff],
           col=col1,
           pch=16,
           cex=0.8)
    points(x=cutoff+1:length(dat),
           y=dat[cutoff+1:length(dat)],
           col=col2,
           pch=16,
           cex=0.8)
}

sing.vals <- read.table("#{ARGV.first}")$V1
var.explained <- sing.vals ^ 2 / sum(sing.vals ^ 2) * 100

pdf("#{var_plot_pdf}", width=8, height=5)
plot.colored.by.inflection.point(var.explained, #{ipoint}, ylab="Variance explained")
invisible(dev.off())
}

tmp_r_fname = File.join File.dirname(ARGV.first), "sv.plot_variance.r"

File.open(tmp_r_fname, "w") do |f|
  f.puts rscript_str
end

Process.run_and_time_it! "Plotting variance", "Rscript #{tmp_r_fname}"

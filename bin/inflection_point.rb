#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require "abort_if"
require "matrix"

include AbortIf

def vector_projection a, b
  (a.dot(b) / b.dot(b).to_f) * b
end

def vector_rejection a, b
  a - vector_projection(a, b)
end

def index_of_max ary
  ary.each_with_index.max.last
end

# Returns the (index + 1) of the inflection point. Doubles as number
# of items until the inflection point.
def inflection_point dat
  if dat.length < 3
    AbortIf.logger.warn { "Fewer than tree data points...." +
                          "Use all of them." }
    dat.length
  else
    num_points = dat.length
    first_point = Vector[0, dat.first]
    last_point = Vector[num_points-1, dat.last]

    b = last_point - first_point

    dists = (0..num_points-1).map do |idx|
      this_point = Vector[idx, dat[idx]]
      a = this_point - first_point
      vector_rejection(a, b).norm
    end

    p dists

    index_of_max(dists) + 1 # add one for the 1.. range
  end
end

data = File.open(ARGV.first, "rt").read.split.map(&:to_f).sort.reverse
p data
AbortIf.logger.info { "Inflection point for sorted, reversed data: #{inflection_point(data)}" }

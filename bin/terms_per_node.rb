#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

require "abort_if"
require "bio"
require "set"
require "trollop"

include AbortIf

def leaf? tree, node
  tree.children(node).empty?
end

opts = Trollop.options do
  banner <<-EOS

  Parses an LSA tree.

  Options:
  EOS

  opt(:tree,
      "The LSA tree to parse",
      type: :string)
end

abort_if opts[:tree].nil?,
         "--tree is a required arg. Try #{__FILE__} --help for help."
abort_unless_file_exists opts[:tree]

str = File.read opts[:tree]
newick = Bio::Newick.new str
tree = newick.tree

num = 0
clusters = []
tree.each_node do |node|
  num += 1

  if leaf? tree, node
    cluster_members = [node]
  else
    cluster_members =
      tree.descendents(node, tree.root).select{|node| leaf?(tree, node) }
  end

  clusters << cluster_members
end

# clusters.sort_by { |cluster| cluster.size }.reverse.each_with_index do |cluster, idx|
#   puts [idx + 1, cluster.size, cluster].join "\t"
# end

clusters.each_with_index do |cluster, idx|
  puts [idx + 1, cluster.size, cluster].join "\t"
end

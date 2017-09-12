#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

MAKE_HEATMAP_FUNC = %Q{
open.png <- function(png.fname, type)
{
  png(png.fname, type = type, width = 8, height = 5, units = "in", res = 300)
}

make.heatmap <- function(topic.idx, rowv=NULL)
{
  topic <- subset(dat, dat$topic == topic.idx)

  counts <- topic[, data.cols]

  log.counts <- log(1 + counts)

  log.counts.mat <- t(as.matrix(log.counts))

  rownames(log.counts.mat) <- doc.names

  ## Order the rows based on the tree. It will have the same names as
  ## the doc.names array. But only if the dendro is valid.
  if (is.null(dendro)) {
    log.counts.mat.ordered <- log.counts.mat
  } else {
    log.counts.mat.ordered <- log.counts.mat[tre$tip.label, ]
  }

  png.fname <- paste(png.fname.base,
                     paste("topic", topic.idx, sep="_"),
                     "png",
                     sep=".")

  ## Try to open the png with various types
  tryCatch({
    open.png(png.fname, "quartz")
  }, error = function(err) {
    tryCatch({
      open.png(png.fname, "cairo-png")
    }, error = function(err) {
      tryCatch({
        open.png(png.fname, "cairo")
      }, error = function(err) {
        ## If this one doesn't work, just let the program die
        open.png(png.fname, "Xlib")
      })
    })
  })

  ## Just use defualt ordering for Rowv if the input dedro is NULL
  if (is.null(rowv)) {
      rowv <- TRUE
  }

  heatmap.2(log.counts.mat.ordered,
            Rowv = rowv,
            trace = "none",
            col=colorRampPalette(c("beige", "cadetblue", "darkblue")),
            cexRow = 0.7,
            cexCol = 0.7,
            labCol = NA, # topic$term.idx,
            adjRow = c(NA, 0.5),
            adjCol = c(NA, 0.5),
            key.title = "",
            key.xlab = "Adjusted count")

  invisible(dev.off())
}
}

def run_rscript! fname
  cmd = "Rscript #{fname}"
  Process.run_and_time_it! "Running RScript", cmd
end

ORIGINAL_MD_GROUP_START_IDX = 3
NEW_MD_GROUP_START_IDX  = 6

# data_col_start_idx is one based (use 3 for the original data, and 6
# for the metadata groups)
def rscript count_mat, doc_tree, data_col_start_idx, png_outbase
  %Q{
library(gplots)
library(ape)

#{MAKE_HEATMAP_FUNC}

tre <- read.tree("#{doc_tree}")

dendro <- tryCatch ({
  as.dendrogram(as.hclust(tre))
}, error = function(err) {
  print(paste("Error with doc dendrogram for (#{png_outbase}), not using it for the heatmap. If there were only two docs, this is not an issue. If there were more than two docs, the row trees on the heatmaps may look different than the doc trees in the metadata groups folders. R Error: ", err))
  NULL
})

png.fname.base <- "#{png_outbase}"


dat <- read.table("#{count_mat}", header=T, sep="\t")
data.cols <- #{data_col_start_idx}:ncol(dat)
doc.names <- names(dat)[data.cols]

topics <- sort(unique(dat$topic))

lapply(topics, make.heatmap, rowv=dendro)

}
end

require "aai"
require "abort_if"
require "fileutils"
require "lsa"


include AbortIf
include Lsa

Process.extend Aai::CoreExtensions::Process

# Command line args
lsa_outdir    = ARGV[0]
mapping_fname = ARGV[1]

def new_term_row total_terms
  Array.new total_terms, 0
end

def mean ary
  ary.reduce(:+) / ary.count.to_f
end

current = ""


original_metadata_dir = File.join lsa_outdir,
                                  "metadata_groups",
                                  "original"
rscript_dir = File.join original_metadata_dir,
                         "r"

td_matrix_dir = File.join original_metadata_dir,
                          "td_matrix"
top_terms_by_topic_dir = File.join original_metadata_dir,
                                   "top_terms_by_topic"

td_matrix_fname = File.join td_matrix_dir,
                            "td_matrix.txt"
idx_to_doc_fname = File.join td_matrix_dir,
                             "td_matrix.idx_to_doc_map.txt"
top_terms_by_topic_fname = File.join top_terms_by_topic_dir,
                                     "top_terms_by_topic.txt"

abort_unless_file_exists td_matrix_fname
abort_unless_file_exists idx_to_doc_fname
abort_unless_file_exists top_terms_by_topic_fname

by_cluster_dir = File.join top_terms_by_topic_dir, "by_cluster"

original_trees_dir = File.join original_metadata_dir, "trees"
original_doc_tree_fname = File.join original_trees_dir, "doc_dist_tree.newick.txt"
abort_unless_file_exists original_doc_tree_fname

if Dir.exist? by_cluster_dir
  AbortIf.logger.warn { "Outdir #{by_cluster_dir} already exists. It's contents will be overwritten." }
end

FileUtils.mkdir_p by_cluster_dir

original_output_fname = File.join by_cluster_dir, "original.txt"
original_output_for_r_fname = File.join by_cluster_dir, "original_for_r.txt"
md_groups_output_fname = File.join by_cluster_dir, "md_groups.txt"

# Make doc to idx map
doc2idx = {}
File.open(idx_to_doc_fname, "rt").each_line do |line|
  idx, doc = line.chomp.split "\t"

  abort_if doc2idx.has_key?(doc),
           "Doc #{doc} is repeated in #{idx_to_doc_fname}"

  doc2idx[doc] = idx.to_i
end
idx2doc = doc2idx.invert


total_terms = File.open(td_matrix_fname, "rt").each_line.count

# First row is the first doc, columns are terms
count_mat = Array.new(1) { new_term_row(total_terms) }

num_rows = 0
num_cols = 0
total_entries = 0
lines = []

md_group_dir = File.join by_cluster_dir, "metadata_groups"
FileUtils.mkdir_p md_group_dir
# TODO make it selectable from mean, or max
md_group_outf_basename = File.join md_group_dir, "means.by_group"
md_group2outf, doc2new_doc =
               parse_mapping_file mapping_fname,
                                  md_group_outf_basename

File.open(td_matrix_fname, "rt").each_line.with_index do |line, term_idx|
  # STDERR.printf "PARSING -- #{ridx}\r" if (ridx % 100).zero?

  line.chomp.split(" ").each do |col|
    doc_idx, val = col.split(":")

    doc_idx = doc_idx.to_i
    val = val.to_f

    num_seen_docs = count_mat.count
    if doc_idx >= num_seen_docs
      count_mat << new_term_row(total_terms)
    end

    count_mat[doc_idx][term_idx] = val
  end
end

total_docs = count_mat.count

topic2top_terms = {}
File.open(top_terms_by_topic_fname, "rt").each_line do |line|
  topic, term_idx, term, weight = line.chomp.split
  topic = topic.to_i
  term_idx = term_idx.to_i
  weight = weight.to_f

  topic2top_terms[topic] = {} unless topic2top_terms.has_key?(topic)
  topic2top_terms[topic][term_idx] = { term: term, weight: weight }
end

term_count_info = {}

term2info = {}

doc_names =
  idx2doc.sort_by { |idx, doc| idx }.map { |idx, doc| doc }

File.open(original_output_for_r_fname, "w") do |f|
  f.puts ["topic", "term.idx", doc_names].join "\t"
end

count_mat_transpose = count_mat.transpose

File.open(original_output_fname, "w") do |f|
  f.puts ["topic",
          "term.idx",
          "doc.idx",
          "count.in.doc",
          "term",
          "md.groups"].join "\t"
  topic2top_terms.each do |topic, top_terms|
    top_terms.each do |term_idx, info|
      term = info[:term]
      weight = info[:weight]

      unless term2info.has_key? term
        term2info[term] = {}
      end

      term2info[term][:weight] = weight
      term2info[term][:topic]  = topic
      term2info[term][:idx] = term_idx

      term_row = count_mat_transpose[term_idx]
      File.open(original_output_for_r_fname, "a") do |rf|
        rf.puts [topic, term_idx, term_row].join "\t"
      end

      total_docs.times do |doc_idx|
        count_in_doc = count_mat[doc_idx][term_idx]


        original_doc = idx2doc[doc_idx]

        doc_name = idx2doc[doc_idx]

        group_info = doc2new_doc[original_doc].map do |group_name, group|
          unless term_count_info.has_key? term
            term_count_info[term] = {}
          end

          unless term_count_info[term].has_key? group_name
            term_count_info[term][group_name] = {}
          end

          unless term_count_info[term][group_name].has_key? group
            term_count_info[term][group_name][group] = []
          end

          term_count_info[term][group_name][group] << count_in_doc

          [group_name, group].join "~"
        end

        f.puts [topic,
                term_idx,
                doc_idx,
                count_in_doc,
                term,
                [["original", doc_name].join("~"),
                 group_info,].join(","),
               ].join "\t"
      end
    end
  end
end

File.open(md_groups_output_fname, "w") do |f|
  f.puts ["term.idx",
          "term",
          "topic",
          "term.weight",
          "group.name",
          "group",
          "mean.count",
          "min.count",
          "max.counts",
          "all.counts"].join "\t"

  # Write the headers
  the_term, the_info = term_count_info.first
  the_info.each do |group_name, group_info|
    md_group2outf[group_name].puts (
      ["term.idx",
       "term",
       "topic",
       "term.weight",
       "abs.term.weight", # for sorting in excel
       group_info.keys,
      ].join "\t"
    )
  end

  term_count_info.each do |term, info1|
    info1.each do |group_name, info2|

      md_group2outf[group_name].puts(
        [
          term2info[term][:idx],
          term,
          term2info[term][:topic],
          term2info[term][:weight],
          term2info[term][:weight].abs,
          info2.map { |group, counts| mean(counts) }
        ].join "\t"
      )

      info2.each do |group, counts|
        f.puts [term2info[term][:idx],
                term,
                term2info[term][:topic],
                term2info[term][:weight],
                group_name,
                group,
                mean(counts),
                counts.min,
                counts.max,
                counts.join(",")].join "\t"
      end
    end
  end
end

md_group2outf.each { |group_name, outf| outf.close }

heatmap_group_dir = File.join md_group_dir, "per_topic_heatmaps"
FileUtils.mkdir_p heatmap_group_dir

md_group2outf.each do |group_name, outf|
  doc_tree_fname = File.join lsa_outdir, "metadata_groups",
                             group_name,
                             "trees",
                             "doc_dist_tree.newick.txt"

  abort_unless_file_exists doc_tree_fname

  png_outbase = File.join heatmap_group_dir,
                          File.basename(outf.path, ".txt")

  rscript_str = rscript File.absolute_path(outf.path),
                        doc_tree_fname,
                        NEW_MD_GROUP_START_IDX,
                        png_outbase

  rscript_fname = File.join rscript_dir,
                            File.basename(outf.path, ".txt") +
                            ".heatmap_script.r"

  File.open(rscript_fname, "w") do |f|
    f.puts rscript_str
  end

  run_rscript! rscript_fname
end


original_heatmaps_rscript_fname =
  File.join rscript_dir,
            File.basename(original_output_for_r_fname, ".txt") +
            ".heatmap_script.r"

png_outbase = File.join heatmap_group_dir,
                        File.basename(original_output_for_r_fname, ".txt")

rscript_str = rscript original_output_for_r_fname,
                      original_doc_tree_fname,
                      ORIGINAL_MD_GROUP_START_IDX,
                      png_outbase

File.open(original_heatmaps_rscript_fname, "w") do |f|
  f.puts rscript_str
end

run_rscript! original_heatmaps_rscript_fname

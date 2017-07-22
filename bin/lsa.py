#!/usr/bin/env python

from gensim import corpora, models, similarities
import logging
import sys
import os
import glob

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)
logger = logging.getLogger()


# Check command line args
if len(sys.argv) != 4:
    sys.exit("\n\nUSAGE: python %s num-topics td-matrix-file idx-to-term-file" % sys.argv[0])

# Command line args
NUM_TOPICS = int(sys.argv[1])
td_matrix_fname = sys.argv[2]
idx2term_fname = sys.argv[3]

# Output file names
transformed_doc_matrix_fname = td_matrix_fname + ".lsa_py.transformed_doc_matrix.txt"
transformed_doc_dist_matrix_fname = td_matrix_fname + ".lsa_py.transformed_doc_dist.txt"
rows_are_terms_fname = td_matrix_fname + ".lsa_py.rows_are_terms.txt"
singular_values_fname = td_matrix_fname + ".lsa_py.singular_values.txt"
top_terms_fname = td_matrix_fname + ".lsa_py.top_terms.txt"
similarity_tmp_fname = td_matrix_fname + ".similarity_tmp"

# Load corpus iterator form a Matrix Market file on disk. Reads one
# doc at a time, with all terms for that doc in memory.
corpus = corpora.MmCorpus(td_matrix_fname)

# Set up idx2term dictionary
idx2term = {}
with open(idx2term_fname, "U") as f:
    for line in f:
        (idx, term) = line.split()
        idx2term[int(idx)] = term

# Train the lsa model on the corpus
lsi = models.LsiModel(corpus=corpus,
                      id2word=idx2term,
                      num_topics=NUM_TOPICS)

## These are the top terms for each of the topics!!!!
topics = lsi.show_topics(num_topics=-1,
                         num_words=20,
                         log=False,
                         formatted=False)

logger.info("Writing top terms")
with open(top_terms_fname, "w") as f:
    f.write("topic term val\n")
    for topic, terms in topics:
        for (term, val) in terms:
            f.write("%d %s %.5f\n" % (topic, term, val))

# These are the singular values arranged from largest to smallest.
logger.info("Writing singular values")
with open(singular_values_fname, "w") as f:
    for val in lsi.projection.s:
        f.write("%.5f\n" % val)

# Each vec in lsi.projection.u contains the weights for each term
# across the topics that are kept. This can be up to num_topics, but
# it may be less if keeping additional topics doesn't add value.

## TODO what we really want are the variable/term loadings...
## i.e. (U %*% S) / sqrt(n - 1)
logger.info("Writing term weights")
with open(rows_are_terms_fname, "w") as f:
    for vec in lsi.projection.u:
        f.write(" ".join([str(round(num, 5)) for num in vec]))
        f.write("\n")

# These are the docs in latent space!!!!
logger.info("Writing docs in topic space")
with open(transformed_doc_matrix_fname, "w") as f:
    idx = -1
    # this is U^-1*X or V*S, when the SVD of X (corpus) is X=U*S*V^T
    # TODO we'd rather have V * sqrt(n - 1) i.e. principal components
    # scaled to unit variance.
    for doc in lsi[corpus]:
        idx += 1
        vals = []
        for (topic, val) in doc:
            vals.append(str(val))

        f.write("%s\n" % " ".join(vals))

# This is the cosine distance matrix for the docs themselves in latent
# space
logger.info("Writing doc in topic space distance matrix")
index = similarities.Similarity(similarity_tmp_fname,
                                lsi[corpus],
                                num_features=NUM_TOPICS)
sims = index[lsi[corpus]]
with open(transformed_doc_dist_matrix_fname, "w") as f:
    for vec in sims:
        # TODO perhaps angular distance is a better way to convert to
        # distance.
        f.write(" ".join([str(round(1 - abs(num), 5)) for num in vec]))
        f.write("\n")

# Remove the similarity tmp files
for fname in glob.glob(similarity_tmp_fname + "*"):
    os.remove(fname)

logger.info("lsa.py done!")

#!/usr/bin/env python

from gensim import corpora, models, similarities
import logging
import sys
import os
import glob

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s',
                    level=logging.INFO)
logger = logging.getLogger()

NUM_TOPICS = 20

td_matrix_fname = sys.argv[1]
transformed_doc_matrix_fname = td_matrix_fname + ".lsa_py.transformed_doc_matrix.txt"
transformed_doc_dist_matrix_fname = td_matrix_fname + ".lsa_py.transformed_doc_dist.txt"
rows_are_terms_fname = td_matrix_fname + ".lsa_py.rows_are_terms.txt"
singular_values_fname = td_matrix_fname + ".lsa_py.singular_values.txt"
top_terms_fname = td_matrix_fname + ".lsa_py.top_terms.txt"
similarity_tmp_fname = td_matrix_fname + ".similarity_tmp"

# Load corpus iterator form a Matrix Market file on disk.
corpus = corpora.MmCorpus(td_matrix_fname)

# print(corpus)

try:
    idx2term = {}

    idx2term_fname = sys.argv[2]
    with open(idx2term_fname, "U") as f:
        for line in f:
            (idx, term) = line.split()
            idx2term[int(idx)] = term

    lsi = models.LsiModel(corpus=corpus,
                          id2word=idx2term,
                          num_topics=NUM_TOPICS)
except:
    lsi = models.LsiModel(corpus=corpus,
                          num_topics=NUM_TOPICS)


print(lsi)

## These are the top terms for each of the topics!!!!
topics = lsi.show_topics(num_topics=-1, num_words=20, log=False, formatted=False)

logger.info("Writing top terms")
with open(top_terms_fname, "w") as f:
    f.write("topic term val\n")
    for topic, terms in topics:
        for (term, val) in terms:
            f.write("%d %s %.4f\n" % (topic, term, val))



# These are the singular values arranged from largest to smallest.
logger.info("Writing singular values")
with open(singular_values_fname, "w") as f:
    for val in lsi.projection.s:
        f.write("%.5f\n" % val)

# Each vec in lsi.projection.u contains the weights for each term
# across the topics that are kept. This can be up to num_topics, but
# it may be less if keeping additional topics doesn't add value.
# idx = 0
logger.info("Writing term weights")
with open(rows_are_terms_fname, "w") as f:
    for vec in lsi.projection.u:
        f.write(" ".join([str(round(num, 5)) for num in vec]))
        f.write("\n")

# These are the docs in latent space!!!!
logger.info("Writing docs in topic space")
with open(transformed_doc_matrix_fname, "w") as f:
    idx = -1
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
        f.write(" ".join([str(round(1 - num, 5)) for num in vec]))
        f.write("\n")

# Remove the similarity tmp files
for fname in glob.glob(similarity_tmp_fname + "*"):
    os.remove(fname)

logger.info("lsa.py done!")

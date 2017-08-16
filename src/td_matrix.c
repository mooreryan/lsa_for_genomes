#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "err_codes.h"
#include "tommy_helper.h"

#include "../vendor/tommyarray.h"
#include "../vendor/tommyhashlin.h"

#define MAX_STR_LEN 10000
#define DOC_SEP "~"

typedef double idf_func_t(int total_docs, int docs_with_term);

/* num_words_in_doc not always the same as number of terms in doc. If
   each term has a count of 1 in that doc then they will be the
   same. But if some terms have a count > 1 in that doc then the
   num_words_in_docs > num_terms_in_doc */
typedef double tf_func_t(int raw_count, int num_words_in_doc);

int
true_singleton(int raw_count, int docs_with_term)
{
  return raw_count == 1 && docs_with_term == 1;
}

/* To drop singletons completely, pass 0 for the weight param. atc
   stands for adjust term count */
int
atc_weight_singletons(int raw_count, int docs_with_term, int weight)
{
  assert(weight >= 0);
  if (true_singleton(raw_count, docs_with_term)) {
    return weight * raw_count;
  } else {
    return raw_count;
  }
}

double
tf(int raw_count, int num_words_in_doc)
{
  return raw_count;
}

double
tf_binary(int raw_count, int num_words_in_doc)
{
  if (raw_count > 0) {
    return 1;
  } else {
    return 0;
  }
}

double
tf_freq(int raw_count, int num_words_in_doc)
{
  assert(num_words_in_doc != 0);
  return raw_count / (double)num_words_in_doc;
}

/* This one should be good for downplaying the importance of crazy
   high count things like tRNAs */
double
tf_log_norm(int raw_count, int num_words_in_doc)
{
  if (raw_count == 0) {
    return 0;
  } else {
    return 1 + log(raw_count);
  }
}

/* TODO write domain assertions for all the idf functions */

double
idf(int total_docs, int docs_with_term)
{
  assert(docs_with_term != 0);
  return log(total_docs / (double)docs_with_term);
}

double
idf_smooth(int total_docs, int docs_with_term)
{
  assert(docs_with_term != 0);
  return log(1 + (total_docs / (double)docs_with_term));
}

/* Docs that are in every term will blow this up. So subtract one */
/* double */
/* idf_prob(int total_docs, int docs_with_term) */
/* { */
/*   int val = 0; */

/*   if (docs_with_term == total_docs) { */
/*     val = total_docs + 1; */
/*   } else { */
/*     val = total_docs; */
/*   } */

/*   assert(val - docs_with_term > 0); */
/*   assert(docs_with_term != 0); */
/*   return log( (val - docs_with_term) / docs_with_term ); */
/* } */

/* Take the two params so that it can be used like the other idf
   functions */
double
idf_const(int total_docs, int docs_with_term)
{
  return 1;
}

double
weight_count(int raw_term_count,
             int num_words_in_doc,
             int total_docs,
             int docs_with_term,
             tf_func_t* tf_func,
             idf_func_t* idf_func)

{
  return tf_func(raw_term_count, num_words_in_doc) *
    idf_func(total_docs, docs_with_term);
}

/* Will exit on memory error. Returns an allocated string. */
char*
join(const char* str1, const char* str2, char sep)
{
  int size =
    strnlen(str1, MAX_STR_LEN) +
    strnlen(str2, MAX_STR_LEN) +
    2; /* one for the null char, one for the sep */

  char* buf =
    malloc(sizeof(*buf) * size);

  int ret = snprintf(buf, size, "%s%c%s", str1, sep, str2);

  /* TODO size or size - 1 */
  PANIC_UNLESS(ret >= 0 && ret < size,
               MEM_ERR,
               stderr,
               "Memory error while joining %s and %s",
               str1,
               str2);

  return buf;
}

/* If str is empty, will return "" in the array */
tommy_array*
tokenize(char* str, char* split_on)
{
  assert(str != NULL);
  assert(split_on != NULL);

  tommy_array* tokens = malloc(sizeof *tokens);
  PANIC_MEM(tokens, stderr);
  tommy_array_init(tokens);

  char* token  = NULL;
  char* string = NULL;
  char* tofree = NULL;
  char* tmp    = NULL;

  tofree = string = strdup(str);
  PANIC_MEM(string, stderr);

  while ((token = strsep(&string, split_on)) != NULL) {
    tmp = strdup(token);
    PANIC_MEM(tmp, stderr);
    tommy_array_insert(tokens, tmp);
  }

  free(tofree);

  return tokens;
}

typedef struct kv_t {
  char* key;
  int val;

  tommy_node node;
} kv_t;

kv_t*
kv_init(char* key, int val)
{
  kv_t* kv = malloc(sizeof *kv);
  PANIC_MEM(kv, stderr);

  kv->key = strdup(key);
  kv->val = val;

  return kv;
}

void
kv_free(kv_t* kv)
{
  free(kv->key);
  free(kv);
}

int
kv_compare(const void* arg, const void* kv)
{
  return strcmp((const char*)arg,
                ((const kv_t*)kv)->key);
}

int
str_eq(const char* s1, const char* s2)
{
  if (strcmp(s1, s2) == 0) {
    return 1;
  } else {
    return 0;
  }
}


/* tommy_foreach_func */
void
kv_zero_out(void* obj)
{
  kv_t* kv = obj;

  kv->val = 0;
}

void
kv_print_map(void* arg, void* data)
{
  kv_t* kv = data;
  FILE* outstream = arg;

  fprintf(outstream,
          "%d\t%s\n",
          kv->val,
          kv->key);
}

int main(int argc, char *argv[])
{
  if (argc != 6) {
    fprintf(stderr,
            "USAGE: %s mmseq2_clusters.tsv outdir tf_func idf_func singleton_weight\n"
            "Note: outdir must already exist.\n",
            argv[0]);

    exit(1);
  }

  char* opt_tf_func  = argv[3];
  char* opt_idf_func = argv[4];
  long opt_singleton_weight = strtol(argv[5], NULL, 10);

  PANIC_UNLESS(opt_singleton_weight >= 0,
               ARG_ERR,
               stderr,
               "singleton weight must be >= 0, got %ld",
               opt_singleton_weight);

  tf_func_t* tf_func = NULL;
  idf_func_t* idf_func = NULL;

  /* Set the tf function */
  if (strncmp(opt_tf_func, "tf", MAX_STR_LEN) == 0) {
    tf_func = tf;
  } else if (strncmp(opt_tf_func, "tf_binary", MAX_STR_LEN) == 0) {
    tf_func = tf_binary;
  } else if (strncmp(opt_tf_func, "tf_freq", MAX_STR_LEN) == 0) {
    tf_func = tf_freq;
  } else if (strncmp(opt_tf_func, "tf_log_norm", MAX_STR_LEN) == 0) {
    tf_func = tf_log_norm;
  } else {
    fprintf(stderr,
            "WARN -- %s is not an option for tf_func."
            " Using tf\n",
            opt_tf_func);

    tf_func = tf;
  }

  /* Set the idf function */
  if (strncmp(opt_idf_func, "idf", MAX_STR_LEN) == 0) {
    idf_func = idf;
  } else if (strncmp(opt_idf_func, "idf_smooth", MAX_STR_LEN) == 0) {
    idf_func = idf_smooth;
  /* } else if (strncmp(opt_idf_func, "idf_prob", MAX_STR_LEN) == 0) { */
  /*   idf_func = idf_prob; */
  } else if (strncmp(opt_idf_func, "idf_const", MAX_STR_LEN) == 0) {
    idf_func = idf_const;
  } else {
    fprintf(stderr,
            "WARN -- %s is not an option for idf_func."
            " Using idf_smooth\n",
            opt_idf_func);

    idf_func = idf_smooth;
  }


  FILE* cluster_tsv = fopen(argv[1], "r");

  PANIC_IF(cluster_tsv == NULL,
           errno,
           stderr,
           "Could not open '%s': %s",
           argv[1],
           strerror(errno));

  double weighted_count = 0.0;

  const char* idx_to_doc_map_suffix  = "td_matrix.idx_to_doc_map.txt";
  const char* idx_to_term_map_suffix = "td_matrix.idx_to_term_map.txt";
  const char* td_matrix_suffix       = "td_matrix.txt";

  char* idx_to_doc_map_fname =
    join(argv[2], idx_to_doc_map_suffix, '/');
  char* idx_to_term_map_fname =
    join(argv[2], idx_to_term_map_suffix, '/');
  char* td_matrix_fname =
    join(argv[2], td_matrix_suffix, '/');

  PANIC_IF_FILE_CAN_BE_READ(stderr, idx_to_doc_map_fname);
  FILE* idx_to_doc_map_f = fopen(idx_to_doc_map_fname, "w");
  PANIC_IF(idx_to_doc_map_f == NULL,
           FILE_ERR,
           stderr,
           "Could not open %s for writing",
           idx_to_doc_map_fname);

  PANIC_IF_FILE_CAN_BE_READ(stderr, idx_to_term_map_fname);
  FILE* idx_to_term_map_f = fopen(idx_to_term_map_fname, "w");
  PANIC_UNLESS(idx_to_term_map_f,
               FILE_ERR,
               stderr,
               "Could not open %s for writing",
               idx_to_term_map_fname);

  PANIC_IF_FILE_CAN_BE_READ(stderr, td_matrix_fname);
  FILE* td_matrix_f = fopen(td_matrix_fname, "w");
  PANIC_UNLESS(td_matrix_f,
               FILE_ERR,
               stderr,
               "Could not open %s for writing",
               td_matrix_fname);

  tommy_array* tokens = NULL;
  char* doc = NULL;

  kv_t* tmp_kv = NULL;

  int doc_idx = 0;
  int term_idx = 0;
  int line_idx = -1;
  int total_lines = 0;

  char* last_term = NULL;

  char* term = malloc(sizeof *term * (MAX_STR_LEN + 1));
  PANIC_MEM(term, stderr);

  char* seq = malloc(sizeof *seq * (MAX_STR_LEN + 1));
  PANIC_MEM(seq, stderr);

  tommy_hashlin* doc2idx = NULL;
  HASHLIN_INIT(doc2idx);

  tommy_hashlin* term2idx = NULL;
  HASHLIN_INIT(term2idx);

  tommy_hashlin* doc_counts = NULL;
  HASHLIN_INIT(doc_counts);

  tommy_hashlin* words_per_doc = NULL;
  HASHLIN_INIT(words_per_doc);

  tommy_array* doc_names = NULL;
  ARRAY_INIT(doc_names);

  char* doc_name = NULL;

  int total_docs = 0;
  int num_terms = 0;

  fprintf(stderr, "LOG -- Starting 1st pass\n");
  /* Read once to get term and doc indices */
  while(fscanf(cluster_tsv, "%s\t%s", term, seq) == 2) {
    ++total_lines;

    tokens = tokenize(seq, DOC_SEP);
    PANIC_UNLESS(tommy_array_size(tokens) == 2,
                 STD_ERR,
                 stderr,
                 "Wrong number of tokens.");

    doc = tommy_array_get(tokens, 0);

    tmp_kv = NULL;
    tmp_kv = tommy_hashlin_search(term2idx,
                                  kv_compare,
                                  term,
                                  tommy_strhash_u32(0, term));

    if (!tmp_kv) { /* new term */
      tmp_kv = kv_init(term, term_idx++);

      tommy_hashlin_insert(term2idx,
                           &tmp_kv->node,
                           tmp_kv,
                           tommy_strhash_u32(0, tmp_kv->key));
    }

    tmp_kv = NULL;
    tmp_kv = tommy_hashlin_search(doc2idx,
                                  kv_compare,
                                  doc,
                                  tommy_strhash_u32(0, doc));

    if (!tmp_kv) { /* new doc */
      tmp_kv = kv_init(doc, doc_idx++);

      tommy_hashlin_insert(doc2idx,
                           &tmp_kv->node,
                           tmp_kv,
                           tommy_strhash_u32(0, tmp_kv->key));

      /* This counting ht will be reset for each term on the second
         read through */
      tmp_kv = kv_init(doc, 0);

      tommy_hashlin_insert(doc_counts,
                           &tmp_kv->node,
                           tmp_kv,
                           tommy_strhash_u32(0, tmp_kv->key));

      tommy_array_insert(doc_names, strdup(doc));
    }

    /* Trak the total number of seqs/words per doc */
    tmp_kv = NULL;
    tmp_kv = tommy_hashlin_search(words_per_doc,
                                  kv_compare,
                                  doc,
                                  tommy_strhash_u32(0, doc));

    if (tmp_kv) {
      ++tmp_kv->val;
    } else {
      tmp_kv = kv_init(doc, 1);
      tommy_hashlin_insert(words_per_doc,
                           &tmp_kv->node,
                           tmp_kv,
                           tommy_strhash_u32(0, tmp_kv->key));
    }

    ARRAY_DONE(tokens, free);
  }
  num_terms = term_idx;
  fclose(cluster_tsv);
  cluster_tsv = fopen(argv[1], "r");
  PANIC_IF(cluster_tsv == NULL,
           errno,
           stderr,
           "Could not open '%s': %s",
           argv[1],
           strerror(errno));

  total_docs = tommy_array_size(doc_names);
  assert(doc_idx == total_docs);

  int docs_with_term = 0;
  int total_non_zero_entries = 0;

  fprintf(stderr, "LOG -- Starting 2nd pass\n");
  /* Read again to get the counts */
  term_idx = 0;
  while(fscanf(cluster_tsv, "%s\t%s", term, seq) == 2) {
    ++line_idx;
    /* Progress */
    if (line_idx % 1000 == 0) {
      fprintf(stderr,
              "Pass 2 -- reading line %d of %d (%.2f%%)\r",
              line_idx,
              total_lines,
              line_idx / (double)total_lines * 100);
    }

    /* First line, so init the last term to the current term */
    if (line_idx == 0) {
      last_term = strdup(term);
    }

    tokens = tokenize(seq, DOC_SEP);
    doc = tommy_array_get(tokens, 0);

    /* TODO BUG this will break if the cluster file is not sorted on
       the first column */
    if (!str_eq(term, last_term)) {
      /* starting a new term */

      /* Get data for idf calc. Scan through each doc to see if this
       * term is present in that doc. */
      for (int i = 0; i < total_docs; ++i) {
        doc_name = tommy_array_get(doc_names, i);

        tmp_kv = NULL;
        tmp_kv = tommy_hashlin_search(doc_counts,
                                      kv_compare,
                                      doc_name,
                                      tommy_strhash_u32(0, doc_name));

        PANIC_UNLESS(tmp_kv,
                     STD_ERR,
                     stderr,
                     "Doc: %s not found in counting hash table",
                     doc);

        /* Doc counts only contains counts for the current term, so if
           it is greater than 0 for a doc, then that term is in the
           doc */
        if (tmp_kv->val > 0) {
          ++docs_with_term;
        }
      }

      int first_non_zero_doc = 1;
      for (int i = 0; i < total_docs; ++i) {
        doc_name = tommy_array_get(doc_names, i);

        tmp_kv = NULL;
        tmp_kv = tommy_hashlin_search(doc_counts,
                                      kv_compare,
                                      doc_name,
                                      tommy_strhash_u32(0, doc_name));

        PANIC_UNLESS(tmp_kv,
                     STD_ERR,
                     stderr,
                     "Doc: %s not found in counting hash table",
                     doc);

        /* TODO will at least one val always be above 0? */
        if (tmp_kv->val > 0) {
          int num_words_in_doc = 0;
          kv_t* word_count = NULL;
          word_count = tommy_hashlin_search(words_per_doc,
                                            kv_compare,
                                            doc_name,
                                            tommy_strhash_u32(0, doc_name));
          PANIC_UNLESS(word_count,
                       STD_ERR,
                       stderr,
                       "Missing doc '%s' from words_per_doc hash table",
                       doc_name);


          num_words_in_doc = word_count->val;
          int adjusted_term_count =
            atc_weight_singletons(tmp_kv->val,
                                  docs_with_term,
                                  opt_singleton_weight);

          weighted_count = weight_count(adjusted_term_count,
                                        num_words_in_doc,
                                        total_docs,
                                        docs_with_term,
                                        tf_func,
                                        idf_func);
          /* this check is to get the spaces right in the output line */
          if (first_non_zero_doc) { /* first doc/col */
            fprintf(td_matrix_f,
                    "%d:%.5f",
                    i,
                    weighted_count);
            first_non_zero_doc = 0;
          } else {
            fprintf(td_matrix_f,
                    " %d:%.5f",
                    i,
                    weighted_count);
          }
        }
      }
      /* finish off the line */
      /* TODO if there were no vals above 0, could this put a blank
         line? */
      fprintf(td_matrix_f, "\n");

      /* Zero out the counting hash */
      tommy_hashlin_foreach(doc_counts,
                            (tommy_foreach_func*)kv_zero_out);

      total_non_zero_entries += docs_with_term;
      docs_with_term = 0;
      ++term_idx;
    } else {
      /* still on the same term */
    }

    /* Make sure the doc for this line is in the doc_counts ht */
    tmp_kv = NULL;
    tmp_kv = tommy_hashlin_search(doc_counts,
                                  kv_compare,
                                  doc,
                                  tommy_strhash_u32(0, doc));

    PANIC_UNLESS(tmp_kv,
                 STD_ERR,
                 stderr,
                 "Doc: %s not found in counting hash table",
                 doc);

    /* This term is found in the current doc, inc the count */
    ++tmp_kv->val;

    ARRAY_DONE(tokens, free);
    free(last_term);
    last_term = strdup(term);
  }
  ++line_idx;
  fprintf(stderr,
          "Pass 2 -- reading line %d of %d (%.2f%%)\n",
          line_idx,
          total_lines,
          line_idx / (double)total_lines * 100);

  /* Catch the last term */
  /* Get data for idf calc */
  for (int i = 0; i < total_docs; ++i) {
    doc_name = tommy_array_get(doc_names, i);

    tmp_kv = NULL;
    tmp_kv = tommy_hashlin_search(doc_counts,
                                  kv_compare,
                                  doc_name,
                                  tommy_strhash_u32(0, doc_name));

    PANIC_UNLESS(tmp_kv,
                 STD_ERR,
                 stderr,
                 "Doc: %s not found in counting hash table",
                 doc);

    if (tmp_kv->val > 0) {
      ++docs_with_term;
    }
  }

  int first_non_zero_doc = 1;
  for (int i = 0; i < total_docs; ++i) {
    doc_name = tommy_array_get(doc_names, i);

    tmp_kv = NULL;
    tmp_kv = tommy_hashlin_search(doc_counts,
                                  kv_compare,
                                  doc_name,
                                  tommy_strhash_u32(0, doc_name));

    PANIC_UNLESS(tmp_kv,
                 STD_ERR,
                 stderr,
                 "Doc: %s not found in counting hash table",
                 doc);

    if (tmp_kv->val > 0) {
      int num_words_in_doc = 0;
      kv_t* word_count = NULL;
      word_count = tommy_hashlin_search(words_per_doc,
                                        kv_compare,
                                        doc_name,
                                        tommy_strhash_u32(0, doc_name));
      PANIC_UNLESS(word_count,
                   STD_ERR,
                   stderr,
                   "Missing doc '%s' from words_per_doc hash table",
                   doc_name);


      num_words_in_doc = word_count->val;
      int adjusted_term_count =
        atc_weight_singletons(tmp_kv->val,
                              docs_with_term,
                              opt_singleton_weight);


      weighted_count = weight_count(adjusted_term_count,
                                    num_words_in_doc,
                                    total_docs,
                                    docs_with_term,
                                    tf_func,
                                    idf_func);

      /* this check is to get the spaces right in the output line */
      if (first_non_zero_doc) {
        fprintf(td_matrix_f,
                "%d:%.5f",
                i,
                weighted_count);
        first_non_zero_doc = 0;
      } else {
        fprintf(td_matrix_f,
                " %d:%.5f",
                i,
                weighted_count);
      }
    }
  }
  /* finish off the line */
  fprintf(td_matrix_f, "\n");

  total_non_zero_entries += docs_with_term;

  fprintf(stderr, "LOG -- writing doc2idx\n");
  /* TODO print out idx + 1 and reverse it */
  tommy_hashlin_foreach_arg(doc2idx, kv_print_map, idx_to_doc_map_f);

  fprintf(stderr, "LOG -- writing term2idx\n");
  tommy_hashlin_foreach_arg(term2idx, kv_print_map, idx_to_term_map_f);

  /* NEXT STEP: print outputs to named files */

  /* Clean up */
  fprintf(stderr, "LOG -- cleaning up\n");
  HASHLIN_DONE(doc_counts, kv_free);
  HASHLIN_DONE(doc2idx, kv_free);
  HASHLIN_DONE(term2idx, kv_free);
  HASHLIN_DONE(words_per_doc, kv_free);

  ARRAY_DONE(doc_names, free);

  free(term);
  free(seq);
  free(last_term);

  free(idx_to_doc_map_fname);
  free(idx_to_term_map_fname);
  free(td_matrix_fname);

  fclose(cluster_tsv);
  fclose(td_matrix_f);
  fclose(idx_to_term_map_f);
  fclose(idx_to_doc_map_f);

  return 0;
}

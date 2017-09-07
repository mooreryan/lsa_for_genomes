#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <zlib.h>

#include "err_codes.h"
#include "tommy_helper.h"

#include "../vendor/kseq.h"
#include "../vendor/tommyarray.h"
#include "../vendor/tommyhashlin.h"

#define READ_BUF_SIZE 4096

KSTREAM_INIT(gzFile, gzread, READ_BUF_SIZE)

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

typedef struct id2info_t {
  char* id;
  char* info;

  tommy_node node;
} id2info_t;

/* does NOT copy the id and char, they should be pre-allocated */
id2info_t*
id2info_init(char* id, char* info)
{
  id2info_t* id2info = malloc(sizeof *id2info);
  PANIC_MEM(id2info, stderr);

  id2info->id = id;
  id2info->info = info;

  return id2info;
}

void
id2info_free(id2info_t* id2info)
{
  free(id2info->id);
  free(id2info->info);
  free(id2info);
}

int
id2info_compare(const void* arg, const void* id2info)
{
  return strcmp((const char*)arg,
                ((const id2info_t*)id2info)->id);
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

#define NUM_MATCHES 2
#define STRING "my apple pie"

char*
regsubstr(char* string, regmatch_t regmatch)
{
  regoff_t substr_len = regmatch.rm_eo - regmatch.rm_so;

  /* add 1 for the null char */
  char* substr = malloc((substr_len + 1) * sizeof(char));
  PANIC_MEM(substr, stderr);

  strncpy(substr,
          &string[regmatch.rm_so],
          substr_len);
  substr[substr_len] = '\0';

  return substr;
}

int main(int argc, char *argv[])
{

  if (argc != 3) {
    fprintf(stderr,
            "USAGE: %s id_to_info.txt search_btab.txt\n",
            argv[0]);

    exit(OPT_ERR);
  }

  char* opt_id_to_info_fname = argv[1];
  char* opt_search_btab_fname = argv[2];

  PANIC_UNLESS_FILE_CAN_BE_READ(stderr, opt_id_to_info_fname);
  PANIC_UNLESS_FILE_CAN_BE_READ(stderr, opt_search_btab_fname);

  tommy_hashlin* ht_id2info = NULL;
  HASHLIN_INIT(ht_id2info);

  gzFile fp;
  kstream_t* kstream;
  kstring_t kstr = {0, 0, 0};
  long line_idx = 0;

  id2info_t* id2info = NULL;

  char* line = NULL;
  char* id = NULL;
  char* desc = NULL;
  size_t size = 0;

  int ret_val = 0;
  regex_t desc_regex;
  regex_t id_regex;
  char regex_err_buf[100];

  regmatch_t subexp[NUM_MATCHES];

  ret_val = regcomp(&desc_regex,
                    "Descriptions=\\[(.*)\\] Members=",
                    REG_EXTENDED);
  if (ret_val != 0) {
    /* regcomp failed */
    regerror(ret_val, &desc_regex, regex_err_buf, sizeof(regex_err_buf));
    fprintf(stderr,
            "ERROR -- error while running desc regcomp: %s\n",
            regex_err_buf);
    exit(ret_val);
  }

  ret_val = regcomp(&id_regex,
                    "^(.*)\t",
                    REG_EXTENDED);
  if (ret_val != 0) {
    /* regcomp failed */
    regerror(ret_val, &id_regex, regex_err_buf, sizeof(regex_err_buf));
    fprintf(stderr,
            "ERROR -- error while running id regcomp: %s\n",
            regex_err_buf);
    exit(ret_val);
  }


  fp = gzopen(opt_id_to_info_fname, "r");
  kstream = ks_init(fp);

  printf("reading: %s\n", opt_id_to_info_fname);
  while (ks_getuntil(kstream, '\n', &kstr, 0) >= 0) {
    if (line_idx++ % 10000 == 0) {
      fprintf(stderr,
              "Reading IDs: %ld\r",
              line_idx);
    }

    /* add the uc prefix if needed */
    if (kstr.s[0] != 'u' && kstr.s[1] != 'c') {
      size = kstr.l + 1 + 2;
      line = malloc(size * sizeof(char));
      PANIC_MEM(line, stderr);
      snprintf(line, size, "uc%s", kstr.s);
    } else {
      size = kstr.l + 1;
      line = malloc(size * sizeof(char));
      PANIC_MEM(line, stderr);
      snprintf(line, size, "%s", kstr.s);
    }


    /* pull out the ID */
    ret_val = regexec(&id_regex, line, NUM_MATCHES, subexp, 0);
    if (ret_val == 0) { /* match */
      id = regsubstr(line, subexp[1]);
    } else if (ret_val == REG_NOMATCH) {
      fprintf(stderr, "NO MATCH\n");
      exit(1);
    } else { /* some other error */
      regerror(ret_val, &id_regex, regex_err_buf, sizeof(regex_err_buf));
      fprintf(stderr, "ERROR -- error while running regexec: %s\n",
              regex_err_buf);
      exit(ret_val);
    }

    /* pull out the descriptions */
    ret_val = regexec(&desc_regex, line, NUM_MATCHES, subexp, 0);
    if (ret_val == 0) { /* match */
      desc = regsubstr(line, subexp[1]);

      /* TODO sometimes the desc ends with a | char, looks bad when
         printed */
    } else if (ret_val == REG_NOMATCH) {
      fprintf(stderr, "NO MATCH\n");
      exit(1);
    } else { /* some other error */
      regerror(ret_val,
               &desc_regex,
               regex_err_buf,
               sizeof(regex_err_buf));
      fprintf(stderr, "ERROR -- error while running regexec: %s\n",
              regex_err_buf);
      exit(ret_val);
    }


    /* abort if the ID is already in the ht */
    id2info = NULL;
    id2info = tommy_hashlin_search(ht_id2info,
                                   id2info_compare,
                                   id,
                                   tommy_strhash_u32(0, id));
    if (id2info) { /* id is repeated in the ht */
      fprintf(stderr,
              "ERROR -- id '%s' is repeated\n",
              id);
      exit(1);
    }

    /* add the descriptions to the ht, keyed on the uc_prefixed ID */
    id2info = id2info_init(id, desc);
    tommy_hashlin_insert(ht_id2info,
                         &id2info->node,
                         id2info,
                         tommy_strhash_u32(0, id2info->id));

    free(line);
  }

  ks_destroy(kstream);
  gzclose(fp);

  fp = gzopen(opt_search_btab_fname, "r");
  kstream = ks_init(fp);

  printf("reading: %s\n", opt_search_btab_fname);
  tommy_array* tokens;
  char* new_id = NULL;

  while (ks_getuntil(kstream, '\n', &kstr, 0) >= 0) {

    tokens = tokenize(kstr.s, "\t");
    PANIC_UNLESS(tommy_array_size(tokens) == 12,
                 STD_ERR,
                 stderr,
                 "Wrong number of tokens for line: %s\n",
                 kstr.s);

    id = tommy_array_get(tokens, 1);
    /* add the uc prefix if needed */
    if (id[0] != 'u' && id[1] != 'c') {
      size = strlen(id) + 1 + 2;
      new_id = malloc(size * sizeof(char));
      PANIC_MEM(new_id, stderr);
      snprintf(new_id, size, "uc%s", id);
    } else {
      size = kstr.l + 1;
      new_id = malloc(size * sizeof(char));
      PANIC_MEM(new_id, stderr);
      snprintf(new_id, size, "%s", id);
    }

    id2info = NULL;
    id2info = tommy_hashlin_search(ht_id2info,
                                   id2info_compare,
                                   new_id,
                                   tommy_strhash_u32(0, new_id));

    /* TODO print out the uc version of the id rather than the
       original */
    if (id2info) { /* the id is present */
      /* print the line + the descriptions */
      fprintf(stdout, "%s\t%s\n", kstr.s, id2info->info);
    } else {
      /* print just the line */
      fprintf(stdout, "%s\n", kstr.s);
    }

    ARRAY_DONE(tokens, free);
    free(new_id);
  }

  ks_destroy(kstream);
  gzclose(fp);

  regfree(&desc_regex);
  regfree(&id_regex);
  free(kstr.s);

  HASHLIN_DONE(ht_id2info, id2info_free);

  return 0;
}

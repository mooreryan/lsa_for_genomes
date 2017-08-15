#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "../vendor/redsvd/vendor/eigen3/Eigen/Dense"
#include "../vendor/redsvd/vendor/eigen3/Eigen/Eigenvalues"

#include "err_codes.h"

/* A guard for snprintf and such */
#define MAX_STR_LEN 10000

/* Will exit on memory error. Returns an allocated string. */
char*
join(const char* str1, const char* str2, char sep)
{
  int size =
    strnlen(str1, MAX_STR_LEN) +
    strnlen(str2, MAX_STR_LEN) +
    2; /* one for the null char, one for the sep */

  char* buf =
    (char*)malloc(sizeof(*buf) * size);
  PANIC_MEM(buf, stderr);

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

namespace PSVD {
  typedef Eigen::Matrix<double,
                        Eigen::Dynamic,
                        Eigen::Dynamic,
                        Eigen::RowMajor> Matrix;
}

/* magnitude of vector */
double
mag(Eigen::VectorXd v)
{
  return sqrt(v.array().square().sum());
}

double
cosine_similarity(Eigen::VectorXd v1, Eigen::VectorXd v2)
{
  return v1.dot(v2) / (mag(v1) * mag(v2));
}

/* TODO other options may be better eg angular distance */
double
cosine_dissimilarity(Eigen::VectorXd v1, Eigen::VectorXd v2)
{
  return 1 - std::abs(cosine_similarity(v1, v2));
}

void
print_all_rows_cos_dis(FILE* fp, PSVD::Matrix &mat)
{
  int num_rows = mat.rows();
  for (int r1 = 0; r1 < num_rows - 1; ++r1) {
    for (int r2 = r1 + 1; r2 < num_rows; ++r2) {
      fprintf(fp,
              "%d %d %+lf\n",
              r1,
              r2,
              cosine_dissimilarity(mat.row(r1), mat.row(r2)));
    }
  }
}


/* Needs to have the first line already read, else it will PANIC */
void
populate_matrix(FILE* fp, int nrows, int ncols, int num_topics, PSVD::Matrix &mat)
{
  int ridx = 0;
  int cidx = 0;
  int ret_val = 0;
  double val = 0.0;

    /* Populate the svd_u matrix */
  for (ridx = 0; ridx < nrows; ++ridx) {
    /* Go until ncols_u to eat all the vals */
    for (cidx = 0; cidx < ncols; ++cidx) {
      ret_val = fscanf(fp, "%lf", &val);
      PANIC_UNLESS(ret_val == 1,
                   STD_ERR,
                   stderr,
                   "Got %d values fp, should have been 1",
                   ret_val);

      /* Don't add all columns, just those requested. */
      if (cidx < num_topics) {
        mat(ridx, cidx) = val;
      }
    }
  }
}

void
print_matrix(FILE* fp, PSVD::Matrix &mat)
{
  for (int i = 0; i < mat.rows(); ++i) {
    fprintf(fp, "%lf", mat(i, 0));

    for (int j = 1; j < mat.cols(); ++j) {
      fprintf(fp, " %lf", mat(i, j));
    }
    fprintf(fp, "\n");
  }
}

int
main(int argc, char *argv[])
{
  if (argc != 6) {
    fprintf(stderr,
            "\nUSAGE: %s svd.U svd.S svd.V num_topics outdir\n"
            "Note: outdir must already exist!\n\n",
            argv[0]);

    return 1;
  }

  int ret_val = 0;

  int nrows_u = 0;
  int ncols_u = 0;
  int nrows_s = 0;
  int ncols_s = 0;
  int nrows_v = 0;
  int ncols_v = 0;

  int ridx = 0;
  int cidx = 0;

  double val = 0.0;

  long num_topics = 0;

  PSVD::Matrix svd_u;
  PSVD::Matrix svd_s;
  PSVD::Matrix svd_v;
  PSVD::Matrix svd_vs;
  PSVD::Matrix svd_us;

  /* Check command line file args */
  PANIC_UNLESS_FILE_CAN_BE_READ(stderr, argv[1]);
  PANIC_UNLESS_FILE_CAN_BE_READ(stderr, argv[2]);
  PANIC_UNLESS_FILE_CAN_BE_READ(stderr, argv[3]);

  FILE* svd_u_f = fopen(argv[1], "r");
  FILE* svd_s_f = fopen(argv[2], "r");
  FILE* svd_v_f = fopen(argv[3], "r");

  PANIC_IF(svd_u_f == NULL,
           errno,
           stderr,
           "Could not open %s for reading. (error: %s)",
           argv[1],
           strerror(errno));
  PANIC_IF(svd_s_f == NULL,
           errno,
           stderr,
           "Could not open %s for reading. (error: %s)",
           argv[2],
           strerror(errno));
  PANIC_IF(svd_v_f == NULL,
           errno,
           stderr,
           "Could not open %s for reading. (error: %s)",
           argv[3],
           strerror(errno));

  num_topics = strtol(argv[4], NULL, 10);
  PANIC_IF(num_topics < 1,
           STD_ERR,
           stderr,
           "Num topics (%ld) must be >= 1",
           num_topics);

  /* Make outfiles */
  char* svd_us_fname = join(argv[5], "svd.US", '/');
  char* svd_vs_fname = join(argv[5], "svd.VS", '/');
  char* svd_vs_dis_fname = join(argv[5], "svd.VS.dis", '/');

  PANIC_IF_FILE_CAN_BE_READ(stderr, svd_us_fname);
  PANIC_IF_FILE_CAN_BE_READ(stderr, svd_vs_fname);
  PANIC_IF_FILE_CAN_BE_READ(stderr, svd_vs_dis_fname);

  FILE* svd_us_f = fopen(svd_us_fname, "w");
  FILE* svd_vs_f = fopen(svd_vs_fname, "w");
  FILE* svd_vs_dis_f = fopen(svd_vs_dis_fname, "w");

  PANIC_IF(svd_us_f == NULL,
           FILE_ERR,
           stderr,
           "Could not open %s for writing",
           svd_us_fname);
  PANIC_IF(svd_vs_f == NULL,
           FILE_ERR,
           stderr,
           "Could not open %s for writing",
           svd_vs_fname);
  PANIC_IF(svd_vs_dis_f == NULL,
           FILE_ERR,
           stderr,
           "Could not open %s for writing",
           svd_vs_dis_fname);

  /* Get num rows and cols for svd_u */
  ret_val = fscanf(svd_u_f, "%d %d", &nrows_u, &ncols_u);

  PANIC_UNLESS(ret_val == 2,
               STD_ERR,
               stderr,
               "Got %d values from the first line of svd_u_f, should be 2",
               ret_val);

  PANIC_UNLESS(num_topics <= ncols_u,
               STD_ERR,
               stderr,
               "num_topics (%ld) must be less than ncols_u (%d)",
               num_topics,
               ncols_u);

  svd_u.resize(nrows_u, num_topics);
  populate_matrix(svd_u_f, nrows_u, ncols_u, num_topics, svd_u);

  /* Read the singular values */
  ret_val = fscanf(svd_s_f, "%d %d", &nrows_s, &ncols_s);

  PANIC_UNLESS(ret_val == 2,
               STD_ERR,
               stderr,
               "Got %d values from the first line of svd_s_f, should be 2",
               ret_val);

  PANIC_UNLESS(nrows_s == ncols_u,
               STD_ERR,
               stderr,
               "nrows_s (%d) did not match ncols_u (%d)",
               nrows_s,
               ncols_u);

  svd_s.resize(num_topics, 1);

  for (ridx = 0; ridx < num_topics; ++ridx) {
    ret_val = fscanf(svd_s_f, "%lf", &val);
    PANIC_UNLESS(ret_val == 1,
                 STD_ERR,
                 stderr,
                 "Got %d values from svd_s_f, should have been 1",
                 ret_val);

    svd_s(ridx, 0) = val;
  }

  /* Read the V matrix */
  ret_val = fscanf(svd_v_f, "%d %d", &nrows_v, &ncols_v);

  PANIC_UNLESS(ret_val == 2,
               STD_ERR,
               stderr,
               "Got %d values from the first line from svd_v_f, should be 2",
               ret_val);

  PANIC_UNLESS(nrows_s == ncols_v,
               STD_ERR,
               stderr,
               "nrows_s (%d) did not match ncols_v (%d)",
               nrows_s,
               ncols_v);

  PANIC_UNLESS(num_topics <= ncols_v,
               STD_ERR,
               stderr,
               "num_topics (%ld) must be less than ncols_v (%d)",
               num_topics,
               ncols_v);

  svd_v.resize(nrows_v, num_topics);
  populate_matrix(svd_v_f, nrows_v, ncols_v, num_topics, svd_v);

  svd_us = svd_u * svd_s.asDiagonal();
  svd_vs = svd_v * svd_s.asDiagonal();

  print_matrix(svd_us_f, svd_us);
  print_matrix(svd_vs_f, svd_vs);

  /* The VS file has the projected docs in our case */
  print_all_rows_cos_dis(svd_vs_dis_f, svd_vs);

  free(svd_us_fname);
  free(svd_vs_fname);
  free(svd_vs_dis_fname);

  fclose(svd_us_f);
  fclose(svd_vs_f);
  fclose(svd_vs_dis_f);

  return 0;
}

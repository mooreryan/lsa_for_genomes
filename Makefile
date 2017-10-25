CC = gcc
MKDIR_P = mkdir -p

CFLAGS = -Wall -g
LDFLAGS = -lz -lm

VENDOR = vendor
SRC = src
BIN = bin
TEST_D = test_files
SMALL_D = $(TEST_D)/small
REAL_D = $(TEST_D)/real
SVD_D = $(VENDOR)/redsvd

OBJS := $(VENDOR)/tommyarray.o \
        $(VENDOR)/tommyhashlin.o \
        $(VENDOR)/tommyhash.o \
        $(VENDOR)/tommylist.o

THREADS = 3

.PHONY: all
.PHONY: clean

all: td_matrix redsvd process_svd

clean:
	-rm -r $(SVD_D)/bin $(SVD_D)/include $(SVD_D)/lib $(SVD_D)/redsvd_build $(BIN)/redsvd
	-rm -r $(BIN)/td_matrix $(BIN)/process_svd $(OBJS)

grep_ids: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

mmseqs_hit_lookup: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

test_mmseqs_hit_lookup: mmseqs_hit_lookup
	$(BIN)/mmseqs_hit_lookup $(TEST_D)/1000_id_to_info.txt $(TEST_D)/test1000.m8

td_matrix: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

test: td_matrix
	rm test_files/td_matrix.*txt; $(BIN)/td_matrix test_files/test.clu.tsv test_files tf idf_smooth 1 && rm test_files/td_matrix.*txt

test_pre:
	rm $(SMALL_D)/*.seanie_lsa.*; ruby prep_seq_files.rb $(THREADS) $(SMALL_D)/*.fa && ruby cluster.rb $(THREADS) $(SMALL_D)/all_clean_annotated.seanie_lsa.fa

test_lsa: td_matrix redsvd process_svd
	rm -r output/; ./lsa.rb -m `which mmseqs` -i test_files/*.faa.gz -a test_files/mapping.txt -t 2 -s 0.95 --grep-seqs

test_docker:
	rm -r output/; bin/run_lsa -i test_files/*.faa.gz -a test_files/mapping.txt -t 2 -s 0.95 --grep-seqs

test_lsa_no_mapping: td_matrix
	rm -r output/; ./lsa.rb -i test_files/*.faa.gz

test_lsa_short: td_matrix process_svd
	rm -r output/metadata_groups; ./lsa.rb -i test_files/*.faa.gz -a test_files/mapping.txt

test_lsa_small:
	rm -r output/; ./lsa.rb -i test_files/small/* -a test_files/mapping.txt

redsvd:
	$(SVD_D)/waf configure --blddir $(SVD_D)/redsvd_build --prefix $(SVD_D)
	$(SVD_D)/waf
	$(SVD_D)/waf install
	mv $(SVD_D)/bin/redsvd ./bin

process_svd:
	g++ src/process_svd.cc -o bin/process_svd

test_process_svd: process_svd
	rm -r test_files/psvd_out; mkdir test_files/psvd_out
	bin/process_svd test_files/svd.U test_files/svd.S test_files/svd.V 2 test_files/psvd_out

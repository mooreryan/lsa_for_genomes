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

MAIN = td_matrix

THREADS = 3

.PHONY: all
.PHONY: clean

all: $(MAIN) redsvd

clean:
	-rm -r $(SVD_D)/bin $(SVD_D)/include $(SVD_D)/lib $(SVD_D)/redsvd_build $(BIN)/redsvd
	-rm -r $(BIN)/td_matrix $(OBJS)

grep_ids: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

test: $(MAIN)
	rm test_files/td_matrix.*txt; $(BIN)/$(MAIN) test_files/test.clu.tsv test_files && rm test_files/td_matrix.*txt

test_pre:
	rm $(SMALL_D)/*.seanie_lsa.*; ruby prep_seq_files.rb $(THREADS) $(SMALL_D)/*.fa && ruby cluster.rb $(THREADS) $(SMALL_D)/all_clean_annotated.seanie_lsa.fa

test_lsa: $(MAIN) redsvd
	rm -r output/; time ./lsa.rb -i test_files/*.faa.gz -a test_files/mapping.txt

test_lsa_no_mapping: $(MAIN)
	rm -r output/; time ./lsa.rb -i test_files/*.faa.gz

test_lsa_short: $(MAIN)
	time ./lsa.rb -i test_files/*.faa.gz -a test_files/mapping.txt

test_lsa_small:
	rm -r output/; time ./lsa.rb -i test_files/small/* -a test_files/mapping.txt

redsvd:
	$(SVD_D)/waf configure --blddir $(SVD_D)/redsvd_build --prefix $(SVD_D)
	$(SVD_D)/waf
	$(SVD_D)/waf install
	mv $(SVD_D)/bin/redsvd ./bin

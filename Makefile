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

OBJS := $(VENDOR)/tommyarray.o \
        $(VENDOR)/tommyhashlin.o \
        $(VENDOR)/tommyhash.o \
        $(VENDOR)/tommylist.o

MAIN = td_matrix

THREADS = 3

.PHONY: all
.PHONY: clean

all: bin_dir $(MAIN) grep_ids

bin_dir:
	$(MKDIR_P) $(BIN)

clean:
	-rm -r $(BIN) $(OBJS) *.o

grep_ids: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

test: $(MAIN)
	rm test_files/test.clu.tsv.*; valgrind --leak-check=full $(BIN)/$(MAIN) test_files/test.clu.tsv && rm test_files/test.clu.tsv.*

test_pre:
	rm $(SMALL_D)/*.seanie_lsa.*; ruby prep_seq_files.rb $(THREADS) $(SMALL_D)/*.fa && ruby cluster.rb $(THREADS) $(SMALL_D)/all_clean_annotated.seanie_lsa.fa

test_lsa: $(MAIN)
	rm -r lsa_output/; time ./lsa.rb -i test_files/*.faa.gz -a test_files/mapping.txt

test_lsa_no_mapping: $(MAIN)
	rm -r lsa_output/; time ./lsa.rb -i test_files/*.faa.gz

test_lsa_short: $(MAIN)
	time ./lsa.rb -i test_files/*.faa.gz -a test_files/mapping.txt

test_lsa_small:
	rm -r lsa_output/; time ./lsa.rb -i test_files/small/* -a test_files/mapping.txt

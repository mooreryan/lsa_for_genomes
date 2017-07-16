CC = gcc
MKDIR_P = mkdir -p

CFLAGS = -Wall -g
LDFLAGS = -lz

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

all: bin_dir $(MAIN)

bin_dir:
	$(MKDIR_P) $(BIN)

clean:
	-rm -r $(BIN) $(OBJS) *.o

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ $(SRC)/$@.c $(LDFLAGS)

test: $(MAIN)
	rm test_files/test.clu.tsv.*; valgrind --leak-check=full $(BIN)/$(MAIN) test_files/test.clu.tsv

test_pre:
	rm $(SMALL_D)/*.seanie_lsa.*; ruby prep_seq_files.rb $(THREADS) $(SMALL_D)/*.fa && ruby cluster.rb $(THREADS) $(SMALL_D)/all_clean_annotated.seanie_lsa.fa

test_pipeline: $(MAIN)
	rm $(REAL_D)/*.seanie_lsa.* $(REAL_D)/*.clean; ruby prep_seq_files.rb $(THREADS) $(REAL_D)/*.fa && ruby cluster.rb $(THREADS) $(REAL_D)/all_clean_annotated.seanie_lsa.fa && $(BIN)/$(MAIN) $(REAL_D)/all_clean_annotated.seanie_lsa.fa.clu.tsv

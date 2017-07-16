# Latent Semantic Analysis For (Meta)Genomes

## Dependencies

The pipleline has the following dependencies.

### RubyGems

- abort_if
- parse_fasta
- trollop
- aai

### Python libraries

- GenSim

### Other programs

- MMseqs2 (expected to be on your `PATH`)

## Installing

Get the code.

```
$ git clone https://github.com/mooreryan/lsa_for_genomes.git
```

Build the C programs.

```
$ cd lsa_for_genomes
$ make
```

You can test that everything works by running

```
$ make test_lsa_full
```

## Running

To see the options, run

```
$ ./lsa.rb -h

  Latent semantic analysis pipeline for genomes and metagenomes.

  --num-topics is the LIMIT on topics. You may see fewer depending on
    the data.

  Options:
  -b, --bin-dir=<s>       Folder with the LSA scripts and binaries (default: ./bin)
  -i, --infiles=<s+>      Files with ORF clusters
  -o, --outdir=<s>        Output directory (default: ./lsa_output)
  -c, --cpus=<i>          Number of CPUs to use (default: 3)
  -n, --num-topics=<i>    The number of topics to calculate for LSA (default: 20)
  -h, --help              Show this message
```

### Example

Here is an example of running the pipeline on the three test files.

```
$ ./lsa.rb -i test_files/*.faa.gz -o lsa_test
```

And you'll get a whole bunch of output...

```
$ tree lsa_test/

lsa_test/
├── all_prepped.fa
├── all_prepped.fa.DB
├── all_prepped.fa.DB.index
├── all_prepped.fa.DB.lookup
├── all_prepped.fa.DB_h
├── all_prepped.fa.DB_h.index
├── all_prepped.fa.clu
├── all_prepped.fa.clu.index
├── all_prepped.fa.clu.tsv
├── all_prepped.fa.clu.tsv.sorted
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.idx_to_doc_map.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.idx_to_term_map.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.rows_are_terms.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.singular_values.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.top_terms.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt.newick.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_matrix.txt
├── lsa_err.txt
├── lsa_out.txt
└── make_tree.r

0 directories, 22 files
```

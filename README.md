# Latent Semantic Analysis For (Meta)Genomes

## Dependencies

The pipleline has the following dependencies.

### Languages

- [Ruby](https://www.ruby-lang.org/es/) (Running the pipeline)
- [Python](https://www.python.org/) (No wonky libraries, just Python itself to run the `waf` script to install `redsvd`)
- [R](https://www.r-project.org/) (For heirarchical clustering)

### RubyGems

- [aai](https://rubygems.org/gems/aai)
- [abort_if](https://rubygems.org/gems/abort_if)
- [lsa](https://rubygems.org/gems/lsa)
- [parse_fasta](https://rubygems.org/gems/parse_fasta)
- [trollop](https://rubygems.org/gems/trollop)

These can be installed from RubyGems like so

```
gem install aii abort_if lsa parse_fasta trollop
```

### R libraries

- [Ape](https://cran.r-project.org/web/packages/ape/index.html)

### Other programs

- [MMseqs2](https://github.com/soedinglab/MMseqs2) (for clustering the ORF files)
- A C and C++ compiler (e.g., [gcc/g++](https://gcc.gnu.org/))

## Installing

Assuming you have met all the dependencies, you need to get the pipeline code. You can clone the repository or download the lastest release from [here](https://github.com/mooreryan/lsa_for_genomes/releases).

```
$ git clone --recursive https://github.com/mooreryan/lsa_for_genomes.git
```

*Note*: don't forget the recursive flag!

Some parts of the pipeline require compiling. You can do that using `make` like so: 

```
$ cd lsa_for_genomes
$ make
```

Test out the pipeline with our toy data set to make sure eveything works okay!

```
$ make test_lsa
```

## Running

To see the options, run

```
$ ./lsa.rb -h

  # Pipeline version: 0.4.0
  # Lib version: 0.1.0
  # Copyright 2017 Ryan Moore
  # Contact: moorer@udel.edu
  # Website: https://github.com/mooreryan/lsa_for_genomes
  # License: MIT

  Latent semantic analysis pipeline for genomes and metagenomes.

  --num-topics is the LIMIT on topics. You may see fewer depending on
    the data.

  Options:
  -b, --bin-dir=<s>             Folder with the LSA scripts and binaries (default: ./bin)
  -m, --mmseqs=<s>              Location of the mmseqs binary (default: /Users/moorer/bin/mmseqs)
  -i, --infiles=<s+>            Files with ORFs
  -o, --outdir=<s>              Output directory (default: ./output)
  -a, --mapping=<s>             Mapping file (optional)
  -c, --cpus=<i>                Number of CPUs to use (default: 4)
  -n, --num-topics=<i>          The maximum number of topics to calculate for LSA (default: 20)
  -p, --percent-of-terms=<i>    What percentage of top terms do you want to look at? (default: 1)
  -v, --version                 Print version and exit
  -h, --help                    Show this message
```

### Example

Here is an example of running the pipeline on the three test files.

```
$ ./lsa.rb -m `which mmseqs` -i test_files/*.faa.gz -o lsa_test
```

*Note*: In this case `mmseqs` was on my path already, so I could do the little backtick trick and pass that to `-m`.

And you'll get a whole bunch of output...

```
$ tree lsa_test/

lsa_output/
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
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.top_terms.txt.ids.txt.topic_1_seqs.fa
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.top_terms.txt.ids.txt.topic_2_seqs.fa
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.top_terms.txt.ids.txt.topic_3_seqs.fa
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_dist.txt.newick.txt
├── all_prepped.fa.clu.tsv.sorted.seanie_lsa.td_matrix.mm.txt.lsa_py.transformed_doc_matrix.txt
├── lsa_err.txt
├── lsa_out.txt
└── make_tree.r

0 directories, 25 files
```

The pipeline is still under active development, so expect this to change.

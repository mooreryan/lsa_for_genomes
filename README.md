# Latent Semantic Analysis For (Meta)Genomes

## Dependencies

The pipleline has the following dependencies. I know...it's a lot :(

### Languages

- [Ruby](https://www.ruby-lang.org/es/)
- [Python](https://www.python.org/)
- [R](https://www.r-project.org/)

### RubyGems

- [abort_if](https://rubygems.org/gems/abort_if)
- [parse_fasta](https://rubygems.org/gems/parse_fasta)
- [trollop](https://rubygems.org/gems/trollop)
- [aai](https://rubygems.org/gems/aai)

These can be installed from RubyGems like so

```
gem install abort_if parse_fasta trollop aai
```

### Python libraries

- [GenSim](https://radimrehurek.com/gensim/)

*Note*: For GenSim to work on one of our servers, we had to upgrade numpy and a Python module called `six`.

### R libraries

- [Ape](https://cran.r-project.org/web/packages/ape/index.html)

### Other programs

- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- A C compiler (e.g., [gcc](https://gcc.gnu.org/))

## Installing

Assuming you have met all the dependencies, you need to get the pipeline code. You can clone the repository or download the lastest release from [here](https://github.com/mooreryan/lsa_for_genomes/releases).

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
$ make test_lsa
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
  -i, --infiles=<s+>      Files with ORFs
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

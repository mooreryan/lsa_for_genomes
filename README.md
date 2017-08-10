# Latent Semantic Analysis For (Meta)Genomes

LSA for Genomes and Metagenomes...it's fun! The pipeline is still under active development so expect things to change quickly. Also, we are working on a preprint with lots of details and usage examples. Stay tuned!

More documentation will be coming soon!

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

If you need to install RubyGems, please go [here](https://rubygems.org/pages/download/).

Using RubyGems, the gems may be installed like so....

```
$ gem install aai abort_if lsa parse_fasta trollop
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

  # Pipeline version: 0.5.0
  # Lib version: 0.1.0
  # Copyright 2017 Ryan Moore
  # Contact: moorer@udel.edu
  # Website: https://github.com/mooreryan/lsa_for_genomes
  # License: MIT

  Latent semantic analysis pipeline for genomes and metagenomes.

  --num-topics is the LIMIT on topics. You may see fewer depending on
    the data. If you pass --num-topics 0, all topics will be used.

  --percent-of-terms-per-topic has an automatic mode, pass 0

  Options:
  -b, --binary-dir=<s>                    Folder with the LSA scripts and binaries (default: /Users/moorer/projects/lsa_for_genomes/bin)
  -m, --mmseqs=<s>                        Location of the mmseqs binary (default: /Users/moorer/bin/mmseqs)
  -i, --infiles=<s+>                      Files with ORFs
  -o, --outdir=<s>                        Output directory (default: /Users/moorer/projects/lsa_for_genomes/output)
  -a, --mapping=<s>                       Mapping file (optional)
  -c, --cpus=<i>                          Number of CPUs to use (default: 4)
  -n, --num-topics=<i>                    The maximum number of topics to calculate for LSA (Use 0 for maximum number of topics.) (default: 0)
  -p, --percent-of-terms-per-topic=<i>    What percentage of top terms per topic do you want to look at? (Use 0 for automatic) (default: 0)
  -v, --version                           Print version and exit
  -h, --help                              Show this message
```

### Example

Here is an example of running the pipeline on the three test files.

```
$ ./lsa.rb -m `which mmseqs` -i test_files/*.faa.gz -o output
```

*Note*: In this case `mmseqs` was on my path already, so I could do the little backtick trick and pass that to `-m`.

And you'll get a whole bunch of output. Here is the structure of the output directory.

```
tree -d -L 3 output/
output/
├── clustering
├── logs
├── metadata_groups
│   ├── metadata_group_1
│   │   ├── r
│   │   ├── redsvd
│   │   ├── td_matrix
│   │   ├── top_terms_by_doc
│   │   ├── top_terms_by_topic
│   │   └── trees
│   ├── metadata_group_2
│   │   ├── r
│   │   ├── redsvd
│   │   ├── td_matrix
│   │   ├── top_terms_by_doc
│   │   ├── top_terms_by_topic
│   │   └── trees
│   └── original
│       ├── r
│       ├── redsvd
│       ├── td_matrix
│       ├── top_terms_by_doc
│       ├── top_terms_by_topic
│       └── trees
└── prepped_seqs

25 directories
```

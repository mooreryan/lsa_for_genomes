# Latent Semantic Analysis For (Meta)Genomes

LSA for Genomes and Metagenomes...it's fun! The pipeline is still under active development so expect things to change quickly. Also, we are working on a preprint with lots of details and usage examples. Stay tuned!

## Installation

### Docker

The LSA pipeline has a fair number of dependencies, so we recommend that you just use the Docker image. You [follow this link](https://www.docker.com/get-docker) to get Docker installed on your computer. See the section on [using LSA](#using-lsa) to see how to pull an image and run LSA with Docker.

### From source

If you are adventerous, you can try to install the LSA pipeline from source. (Again, we *strongly* suggest just using the Docker image).

#### Dependencies

The pipleline has the following dependencies.

##### Languages

You'll need to have interpreters for the following languages installed on your machine.

- [Ruby](https://www.ruby-lang.org/es/) (Running the pipeline)
- [Python](https://www.python.org/) (No wonky libraries, just Python itself to run the `waf` script to install `redsvd`. *Note*: We have found MiniConda not to work with `waf`, but the system Python works fine.)
- [R](https://www.r-project.org/) (For lots of random stuff)

##### RubyGems

- [aai](https://rubygems.org/gems/aai)
- [abort_if](https://rubygems.org/gems/abort_if)
- [bio](https://rubygems.org/gems/bio)
- [bundler](https://rubygems.org/gems/bundler)
- [lsa](https://rubygems.org/gems/lsa)
- [parse_fasta](https://rubygems.org/gems/parse_fasta)
- [trollop](https://rubygems.org/gems/trollop)

If you need to install RubyGems, please go [here](https://rubygems.org/pages/download/).

These will be installed with Bundler. See [the section on installation](#installing).

##### R libraries

- [Ape](https://cran.r-project.org/web/packages/ape/index.html)
- [gplots](https://cran.r-project.org/web/packages/gplots/index.html)

See the documentation for how to install them, but you likely can run something like this

```
$ Rscript -e "install.packages('gplots', repos = 'http://cran.us.r-project.org')"
$ Rscript -e "install.packages('ape', repos = 'http://cran.us.r-project.org')"
```

##### Other programs

- A C and C++ compiler (e.g., [gcc/g++](https://gcc.gnu.org/))
- [GNU Make](https://www.gnu.org/software/make/)
- [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- [MMseqs2](https://github.com/soedinglab/MMseqs2) (for clustering the ORFs)

#### Get the pipeline code

Assuming you have met all the dependencies, you need to get the pipeline code.

##### Using git

You can clone the repository like this:

```
$ git clone --recursive https://github.com/mooreryan/lsa_for_genomes.git
```

*Note*: don't forget the recursive flag!

This will give you the master branch which is continuously updated. If you want a stable release, things are a bit more complicated since the git repo has submodules.

##### Downloading a stable release

First, go and [download the lastest release](https://github.com/mooreryan/lsa_for_genomes/releases).

Then untar the file.

```
$ tar xzf lsa_for_genomes-0.11.4.tar.gz
```

Next, you need to get the code for calculating the SVD, because the release doesn't contain the `redsvd` submodule.

```
$ cd lsa_for_genomes-0.11.4/vendor
$ git clone https://github.com/mooreryan/redsvd.git redsvd/
```

#### Preparing the code to run

##### Install the ruby gems

LSA needs specific versions of certain ruby gems. First you must have the bundler gem installed. Install it like so:

```
$ gem install bundler
```

Then in the source directory run

```
$ bundle install
```

which will manage the ruby gems for you.

*Note*: If you have other versions of the required ruby gems installed, this may break other ruby programs you have. In this case, I recommend you use something like [RVM](https://rvm.io) to manage Ruby and various gemsets.

##### Compiling C and C++ code

Some parts of the pipeline require compiling. You can do that using `make`. From the source directory run

```
$ make
```

##### Test it out

Test out the pipeline with our toy data set to make sure eveything works okay!

```
$ make test_lsa
```

If everything goes well, this will output a folder called `output` in the source directory. Check it out!

## Using LSA

### Prepping your data

The first thing you need to do is to call ORFs on your genomes or contigs. You should have a file of ORFs for the most granular level of analysis that you want to do. For example, one file of ORFs per genome, one file of ORFs per sample, or one file of ORFs per contig.

### Mapping files

Then you can make a metadata mapping file to do higher level groupings of data. The mapping file is a tab delimited text file with a header line. The first column must be "file name" and match the file names of your input files (without the directory part, e.g., if your file is `/home/mooreryan/orfs.faa`, then only put the `orfs.faa` part in column 1.) You can have as many additional columns as you wish.

You can see an example [here](https://raw.githubusercontent.com/mooreryan/lsa_for_genomes/master/test_files/mapping.txt). This mapping file is for the three test genomes.

### Running the pipeline

#### Docker

If you have the docker version, you can use our little helper script to make things easy to run. Here is an example of running the test genomes we provide.

```
bin/run_lsa -i test_files/*.faa.gz -o output
```

This will pull the latest Docker image so that you're up to date, and then run the pipeline using the Docker image. If you don't want to update your Docker image, please use the `docker run` command as described in the Docker tutorials.

#### From source

If you installed from source, here is an example of running the pipeline on the three test files.

```
$ ./lsa.rb -m `which mmseqs` -i test_files/*.faa.gz -o output
```

*Note*: In this case `mmseqs` was on my path already, so I could do the little backtick trick and pass that to `-m`.

And you'll get a whole bunch of output.

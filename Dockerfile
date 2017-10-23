FROM ruby:2.4.2-stretch
MAINTAINER Ryan Moore <moorer@udel.edu>

RUN apt-get -y update && apt-get -y install r-base libopenblas-base gcc g++ cmake git vim

# Install R libs
RUN Rscript -e "install.packages('gplots', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('ape', repos = 'http://cran.us.r-project.org')"

# Install MMseqs2
WORKDIR /opt/mmseqs

# Copying the local mmseqs so we don't have to worry about a new version breaking LSA
COPY vendor/MMseqs2 .

WORKDIR build_sse
RUN cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN make && make install

WORKDIR ../build_avx
RUN cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN make && make install

# Copy mmseqs
RUN mv /opt/mmseqs/build_sse/bin/mmseqs /usr/local/bin/mmseqs_sse42
RUN mv /opt/mmseqs/build_avx/bin/mmseqs /usr/local/bin/mmseqs_avx2
RUN mv /opt/mmseqs/util/mmseqs_wrapper.sh /usr/local/bin/mmseqs

# Set up Ruby and the lsa program

RUN gem install bundler

RUN \curl -sSL https://github.com/mooreryan/lsa_for_genomes/archive/v0.11.3.tar.gz \
    | tar -v -C /home -xz
RUN mv /home/lsa_for_genomes-0.11.3 /home/lsa_for_genomes

WORKDIR /home/lsa_for_genomes
RUN bundle install

WORKDIR /home/lsa_for_genomes/vendor
RUN git clone https://github.com/mooreryan/redsvd.git

WORKDIR /home/lsa_for_genomes
RUN make clean && make

CMD ["ruby", "/home/lsa_for_genomes/lsa.rb", "--help"]

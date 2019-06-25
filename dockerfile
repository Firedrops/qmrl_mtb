FROM google/cloud-sdk
LABEL maintainer="Larry Cai <larrycai.jpl@gmail.com>"

EXPOSE 80 443 22 9418
#22 for SSH, 9418 for git, 80 and 443 general catch-alls

#Uncomment if deployment scenario does not have access to port 9418 which prevents git clones. This is a workaround.
#RUN git config --global url.https://github.com/.insteadOf git://github.com/

RUN apt-get -y update
ENV IMAGE_PACKAGES="zlib1g-dev libz-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev gsalib build-essential make g++ gcc perl python3-pip ggplot2 reshape gplots autoconf automake jq ruby apache2 bwa gzip kalign tar wget vim bedtools r-base libopenblas-base"
RUN apt-get -y install $IMAGE_PACKAGES

#Install Perl CGI module, it's not included into the standard distribution anymore
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install CGI

#Install Samtools
RUN git clone https://github.com/samtools/htslib.git
RUN cd htslib && make
RUN cd /
RUN git clone git://github.com/samtools/samtools.git
RUN cd samtools
RUN autoheader && autoconf -Wno-syntax && ./configure
RUN make
RUN make install
RUN apt-get install pip
RUN cd /

#Install bcftools
RUN git clone git://github.com/samtools/bcftools.git
RUN cd bcftools
RUN autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
RUN make
RUN cd /

#Install Picard
RUN git clone https://github.com/broadinstitute/picard.git
RUN cd picard/
RUN ./gradlew shadowJar
RUN ./gradlew clean
RUN cd /

#Install GATK4
RUN git clone https://github.com/broadinstitute/gatk.git
RUN cd gatk/
RUN ./gradlew bundle
RUN ./gradlew clean
RUN cd /

#Install Trimmomatic
RUN git clone https://github.com/timflutre/trimmomatic.git
RUN cd trimmomatic/
RUN make
RUN check
RUN make install
RUN cd /

#Install FASTQC
RUN git clone https://github.com/s-andrews/FastQC.git
RUN cd FastQC/
RUN chmod 755 fastqc
RUN sudo ln -s /FastQC/fastqc /usr/local/bin/fastqc
RUN cd /

#Install SPAdes // not very elegant, posted issue on github, check back soon
RUN wget http://cab.spbu.ru/files/release3.13.1/SPAdes-3.13.1-Linux.tar.gz
RUN tar -xzf SPAdes-3.13.1-Linux.tar.gz
RUN rm SPAdes-3.13.1-Linux.tar.gz
RUN mv SPAdes-* SPAdes
RUN PATH=$PATH:/SPAdes/bin/
RUN cd /

#Install Freebayes
RUN git clone --recursive git://github.com/ekg/freebayes.git
RUN cd freebayes/
RUN make
RUN make install
RUN cd /

#Install yaggo (prerequisite for Mummer)
RUN gem install yaggo

#Install Mummer NOT DONE YET CHECK GITHUB ISSUES https://github.com/mummer4/mummer/issues/107
RUN git clone https://github.com/mummer4/mummer.git
RUN cd mummer/
RUN autoreconf -fi
RUN ./configure --prefix=/
RUN make
RUN make install
RUN cd /

#Install Kraken2
RUN git clone https://github.com/DerrickWood/kraken2.git
RUN cd kraken2/
RUN ./install_kraken2.sh /kraken2/
RUN cp /kraken2/kraken2{,-build,-inspect} /bin
RUN cd /

#Install beast 1.x
RUN curl -s "https://api.github.com/repos/beast-dev/beast-mcmc/releases/latest" | jq --arg PLATFORM_ARCH "tgz" -r '.assets[] | select(.name | endswith($PLATFORM_ARCH)).browser_download_url' | xargs curl -L -o /beast1.tgz
RUN tar -zxvf beast1.tgz
RUN mv BEAST* beast1
RUN rm beast1.tgz
RUN export PATH=$Path:/beast1/bin/
RUN cd /

#Install beast 2.x #NOT DONE YET CHECK GITHUB ISSUES https://gist.github.com/lukechilds/a83e1d7127b78fef38c2914c4ececc3c
RUN curl -s "https://api.github.com/repos/CompEvol/beast2/releases/latest" | jq --arg PLATFORM_ARCH "BEAST.v2*tgz" -r '.assets[] | select(.name | endswith($PLATFORM_ARCH)).browser_download_url' | xargs curl -L -o /beast2.tgz
RUN tar -zxvf beast1.tgz
RUN mv BEAST* beast1
RUN rm beast1.tgz
RUN export PATH=$Path:/beast1/bin/
RUN cd /

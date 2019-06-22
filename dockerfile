FROM google/cloud-sdk
LABEL maintainer="Larry Cai <larrycai.jpl@gmail.com>"

EXPOSE 80 443 22 9418
#22 for SSH, 9418 for git, 80 and 443 general catch-alls

#Uncomment if deployment scenario does not have access to port 9418 which prevents git clones. This is a workaround.
#RUN git config --global url.https://github.com/.insteadOf git://github.com/

RUN apt-get -y update
ENV IMAGE_PACKAGES="zlib1g-dev libz-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev gsalib build-essential make gcc perl python3-pip ggplot2 reshape gplots autoconf automake apache2 bwa gzip kalign tar wget vim bedtools r-base libopenblas-base"
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
RUN ./gradlew bundle
RUN ./gradlew clean
RUN cd /

#Install Trimmomatic
RUN git clone https://github.com/timflutre/trimmomatic.git
RUN make
RUN check
RUN make install
RUN cd /

#Install FASTQC
RUN git clone https://github.com/s-andrews/FastQC.git
RUN chmod 755 fastqc
RUN sudo ln -s /FastQC/fastqc /usr/local/bin/fastqc
RUN cd /

#Install SPAdes // not very elegant, posted issue on github, check back soon
RUN mkdir -p /SPAdes
RUN wget http://cab.spbu.ru/files/release3.13.1/SPAdes-3.13.1-Linux.tar.gz
RUN tar -xzf SPAdes-3.13.1-Linux.tar.gz
RUN rm SPAdes-3.13.1-Linux.tar.gz
RUN PATH=$PATH:~/SPAdes-3.13.1-Linux/bin/
RUN cd /

#Install Freebayes
RUN git clone --recursive git://github.com/ekg/freebayes.git
RUN make
RUN make Install
RUN cd /

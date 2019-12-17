FROM continuumio/miniconda3
LABEL maintainer="Larry Cai <larrycai.jpl@gmail.com>"

#Install various dependencies. I am not sure that all of them are needed, just cloned all for now to save time.
RUN apt-get -qqy update && apt-get install -qqy \
        curl \
        gcc \
        python-dev \
        python-setuptools \
        apt-transport-https \
        lsb-release \
        openssh-client \
        git \
        gnupg

RUN mkdir -p /usr/share/man/man1
ENV IMAGE_PACKAGES="apt-utils zlib1g-dev libz-dev libbz2-dev liblzma-dev libperl-dev libcurl4-gnutls-dev libgconf-2-4 libssl-dev libncurses5-dev libopenblas-base libtool libx11-dev libxext-dev libxrender-dev libxrandr-dev libxtst-dev libxt-dev libcups2-dev libasound2-dev pkg-config build-essential perl python-pip python3-pip autoconf automake jq ruby apache2 bwa gzip kalign tar wget vim ant bedtools"
RUN apt -y update && apt -y install $IMAGE_PACKAGES
ENV PIP_PACKAGES="gsalib reshape"
RUN pip install $PIP_PACKAGES

RUN curl -L https://cpanmin.us | perl - App::cpanminus \
  && cpanm CGI \
    Statistics::Basic \
    MCE

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt install -y git-lfs

#Install mtbseq(+samtools, picard, bwa)
RUN conda install -c bioconda mtbseq

#Download GATK3, prerequisite for mtbseq. Yes, it should be done AFTER mtbseq install.
RUN mkdir /gatk3
RUN wget -O /gatk3/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 -U "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36" --referer https://software.broadinstitute.org/gatk/download/archive 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef' && \
    gatk3-register /gatk3/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

#Install VCFtools
RUN git clone https://github.com/vcftools/vcftools.git && \
	cd vcftools && \
	./autogen.sh && \
	./configure && \
	make && \
	make install && \
	cd /

#Install beast 2.x
RUN curl -s "https://api.github.com/repos/CompEvol/beast2/releases" | grep download | grep tgz | head -n 1 | awk '{print $2}' | xargs curl -L -o /beast2.tgz
RUN tar -zxvf beast2.tgz && rm beast2.tgz && mv beast beast2
ENV PATH $PATH:/beast2/bin

#Install GSL & htslib, prerequisites for bcftools
RUN git clone git://git.savannah.gnu.org/gsl.git && \
	cd /gsl && \
	./autogen.sh && \
	./configure && \
	make && \
	make install && \
	cd /

RUN git clone https://github.com/samtools/htslib.git && cd htslib && make && make install && cd /

#Install bcftools
RUN git clone git://github.com/samtools/bcftools.git && cd bcftools && make && make install && cd /

#Install Trimmomatic
RUN git clone https://github.com/timflutre/trimmomatic.git && cd trimmomatic/ && make && make check && make install && cd /

#Install SnpEff v4.1l (SnpSift included) outdated but for compatibility
#RUN wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
RUN curl -L https://api.github.com/repos/pcingola/SnpEff/zipball/v4.1l -o /snpEff-v4.1l.zip
RUN unzip /snpEff-v4.1l.zip && rm /snpEff-v4.1l.zip && mv /pcingola* /SnpEff
ENV PATH $PATH:/SnpEff
#Optional: The database directory can be changed in snpEff.config. Default is in the installation folder (./data/).
#See: http://snpeff.sourceforge.net for more information

#Set PERL thingy for vcftools
ENV PERL5LIB=/vcftools/src/perl/

#Make data folder to mount nectar data to
RUN mkdir /data

#Copy scripts
COPY pipe.sh pipe2.sh vcf_filter_module.py vcf_tab_to_fasta_alignment.pl ./

#Make executable
RUN chmod +x /pipe.sh /pipe2.sh /vcf_filter_module.py /vcf_tab_to_fasta_alignment.pl

#Report Version
RUN echo "Version 2019_17_12, revert snpEff to version 4.3l"

WORKDIR /
ENTRYPOINT ["./pipe.sh"]

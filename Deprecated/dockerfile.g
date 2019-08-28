#Based on google cloud docker for built-in integration, does not contain R.

FROM google/cloud-sdk
LABEL maintainer="Larry Cai <larrycai.jpl@gmail.com>"

EXPOSE 80 443 22 9418
#22 for SSH, 9418 for git, 80 and 443 general catch-alls

#Uncomment if deployment scenario does not have access to port 9418 which prevents git clones. This is a workaround.
#RUN git config --global url.https://github.com/.insteadOf git://github.com/

#Install apt packages, mostly dependencies, some are software
RUN apt -y update
ENV IMAGE_PACKAGES="zlib1g-dev libz-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libgconf-2-4 libssl-dev libncurses5-dev libopenblas-base build-essential make g++ gcc perl python3-pip autoconf automake r-base jq ruby apache2 bwa gzip kalign tar wget vim bedtools"
RUN apt -y install $IMAGE_PACKAGES
ENV PIP_PACKAGES="gsalib reshape"
RUN pip install $PIP_PACKAGES

#Install Perl CGI module, it's not included into the standard distribution anymore
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install CGI

#Install other Perl modules
RUN cpanm Statistics::Basic
RUN cpanm MCE

#Install git lfs, prerequisite for gradlew installations
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt install -y git-lfs

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
RUN cd /ls

#Install IGV Possible Debug. My test VM crashes at ./gradlew test but no errors thrown, so assumed working.
RUN git clone https://github.com/igvteam/igv.git
RUN cd igv/
RUN ./gradlew createDist
RUN ./gradlew createToolsDist
RUN ./gradlew test --no-daemon
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

#Install SPAdes
RUN curl -s "https://api.github.com/repos/ablab/spades/releases" | grep download | grep Linux.tar.gz | head -n 1 | awk '{print $2}' | xargs curl -L -o /SPAdes.tar.gz
RUN tar -zxvf SPAdes.tar.gz
RUN rm SPAdes.tar.gz
RUN mv SPAdes* SPAdes
RUN export PATH=$PATH:/SPAdes/bin/
RUN cd /

#Install Freebayes
RUN git clone --recursive git://github.com/ekg/freebayes.git
RUN cd freebayes/
RUN make
RUN make install
RUN cd /

#Install yaggo (prerequisite for Mummer)
RUN gem install yaggo

#Install Mummer
RUN curl -s "https://api.github.com/repos/mummer4/mummer/releases" | grep download | grep tar.gz | head -n 1 | awk '{print $2}' | xargs curl -L -o /mummer.tar.gz
RUN tar -zxvf mummer.tar.gz
RUN rm mummer.tar.gz
RUN mv mummer* mummer
RUN cd mummer/
RUN ./configure --prefix=/mummer/
RUN make
RUN make install
RUN export PATH=$PATH:/mummer/
RUN cd /

#Install Prodigal (Prerequisite for Circlator)
RUN git clone https://github.com/hyattpd/Prodigal.git
RUN cd Prodigal
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
#RUN mv BEAST* beast1
RUN rm beast1.tgz
RUN export PATH=$PATH:/beast1/bin/
RUN cd /

#Install beast 2.x
RUN curl -s "https://api.github.com/repos/CompEvol/beast2/releases/latest" | grep download | grep tgz | head -n 1 | awk '{print $2}' | xargs curl -L -o /beast2.tgz
RUN tar -zxvf beast2.tgz
#RUN mv beast beast2
RUN rm beast2.tgz
RUN export PATH=$PATH:/beast2/bin/
RUN cd /

#Install Figtree
RUN curl -s "https://api.github.com/repos/rambaut/figtree/releases/latest" | jq --arg PLATFORM_ARCH "tgz" -r '.assets[] | select(.name | endswith($PLATFORM_ARCH)).browser_download_url' | xargs curl -L -o /figtree.tgz
RUN tar -zxvf figtree.tgz
RUN mv FigTree* FigTree
RUN rm figtree.tgz
RUN export PATH=$PATH:/FigTree/bin/
RUN cd /

#Install SnpEff (SnpSift included) #To be tested
RUN wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip
RUN rm snpEff_latest_core.gzip
RUN mv snpEff* snpEff
#Optional: The database directory can be changed in snpEff.config. Default is in the installation folder (./data/).
#See: http://snpeff.sourceforge.net for more information

#Install Circos
RUN wget -O /Circos.tgz http://circos.ca/distribution/circos-current.tgz
RUN tar -zxvf Circos.tgz
RUN rm Circos.tgz

#Install Miniconda (Prerequisite for MTBseq)
RUN wget -O /miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash /miniconda.sh -b -f -p /miniconda/
RUN rm /miniconda.sh
RUN export PATH=$PATH:/miniconda/bin/
RUN conda install anaconda
RUN cd /

#Install MTBseq #Almost done check https://github.com/ngs-fzb/MTBseq_source/issues/29 #Still missing depndencies (cpanm stuff)
RUN conda install -y -c bioconda mtbseq
RUN cd /
RUN mkdir /miniconda/dependencies/
RUN wget -O /miniconda/dependencies/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef -U "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36" chromium --referer https://software.broadinstitute.org/gatk/download/archive 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef'
wget -U "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36" chromium --referer https://software.broadinstitute.org/gatk/download/archive 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef'
RUN gatk3-register /miniconda/dependencies/GenomeAnalysisTK[-$PKG_VERSION.tar.bz2|.jar]

#Install Mykrobe #Almost done check https://github.com/Mykrobe-tools/mykrobe/issues/62
RUN conda install -y -c bioconda mykrobe

#Install TBProfiler
RUN conda install -y -c bioconda tb-profiler

#Install Pilon
RUN conda install -y -c bioconda pilon

#Install VCFtools
RUN git clone https://github.com/vcftools/vcftools.git
RUN cd vcftools
RUN ./autogen.sh
RUN ./configure
RUN make
RUN make install
RUN cd /

#Install trimAl
RUN git clone https://github.com/scapella/trimal.git
RUN cd trimal/source
RUN make
RUN export PATH=$PATH:/trimal/source/
RUN cd /

#Install Racon
RUN git clone --recursive https://github.com/isovic/racon.git racon
RUN cd racon
RUN mkdir build
RUN cd build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make
RUN make install
RUN cd /

#Install Canu
RUN curl -s "https://api.github.com/repos/marbl/canu/releases/latest" | jq --arg PLATFORM_ARCH "Linux-amd64.tar.xz" -r '.assets[] | select(.name | endswith($PLATFORM_ARCH)).browser_download_url' | xargs curl -L -o /canu.tar.xz
RUN tar -zJf canu.tar.xz
RUN mv canu-* canu
RUN rm canu.tar.xz
RUN export PATH=$PATH:/canu/Linux-amd64/bin
RUN cd /

#Install Circlator
RUN pip3 install circlator && cd /

#Install QUAST
RUN conda install -y -c bioconda quast && cd /

#Install TempEst
RUN wget -O /tempest.tgz 'http://tree.bio.ed.ac.uk/download.php?id=102&num=3'
RUN tar -zxvf tempest.tgz && rm tempest.tgz
RUN mv TempEst* TempEst
RUN export PATH=$PATH:/TempEst/bin/
RUN cd /

#Install MEGA5 #May become outdated
RUN wget https://www.megasoftware.net/do_force_download/megax_10.0.5-1_amd64.deb #GUI version
RUN dpkg -i megax_10.0.5-1_amd64.deb
RUN cd /

#FOR R: install.packages("tidyverse") includes ggplot2
#install.packages("reshape2")
#install.packages("gplots")
#install.packages("zoo")
#install.packages("cluster")

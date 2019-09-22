#2-in-1 version, default.
#Based on google cloud docker for built-in integration, R installed at the end. See "rstudio series" for notes.

FROM rocker/verse
LABEL maintainer="Larry Cai <larrycai.jpl@gmail.com>"


#Install Google Cloud SDK https://hub.docker.com/r/google/cloud-sdk/dockerfile
ARG CLOUD_SDK_VERSION=263.0.0
ENV CLOUD_SDK_VERSION=$CLOUD_SDK_VERSION
ENV PATH "$PATH:/opt/google-cloud-sdk/bin/"
#COPY --from=static-docker-source /usr/local/bin/docker /usr/local/bin/docker
RUN apt-get -qqy update && apt-get install -qqy \
        curl \
        gcc \
        python-dev \
        python-setuptools \
        apt-transport-https \
        lsb-release \
        openssh-client \
        git \
        gnupg \
    && easy_install -U pip && \
    pip install -U crcmod   && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && \
    apt-get install -y google-cloud-sdk=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-app-engine-python=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-app-engine-python-extras=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-app-engine-java=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-app-engine-go=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-datalab=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-datastore-emulator=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-pubsub-emulator=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-bigtable-emulator=${CLOUD_SDK_VERSION}-0 \
        google-cloud-sdk-cbt=${CLOUD_SDK_VERSION}-0 \
        kubectl
	
VOLUME ["/root/.config"]


EXPOSE 80 443 22 9418
#22 for SSH, 9418 for git, 80 and 443 general catch-alls

#Uncomment if deployment scenario does not have access to port 9418 which prevents git clones. This is a workaround.
#RUN git config --global url.https://github.com/.insteadOf git://github.com/

#Install apt packages, mostly dependencies, some are software
ENV IMAGE_PACKAGES="zlib1g-dev libz-dev libbz2-dev liblzma-dev libperl-dev libcurl4-gnutls-dev libgconf-2-4 libssl-dev libncurses5-dev libopenblas-base libtool pkg-config build-essential perl python-pip python3-pip autoconf automake jq ruby apache2 bwa gzip kalign tar wget vim ant bedtools"
RUN apt -y update && apt -y install $IMAGE_PACKAGES
ENV PIP_PACKAGES="gsalib reshape"
RUN pip install $PIP_PACKAGES

#Install Perl CGI, MCE, and Statistics::Basic modules.
RUN curl -L https://cpanmin.us | perl - App::cpanminus \
  && cpanm CGI \
    Statistics::Basic \
    MCE

#Install git lfs, prerequisite for gradlew installations
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt install -y git-lfs

#Install GSL, prerequisite for bcftools
#RUN git clone https://github.com/LuaDist/gsl.git && \
#	cd /gsl && \
#	./configure --disable-shared && \
#	make && \
#	make install && \
#	cd /

#Install GSL, prerequisite for bcftools ALTERNATIVE	
RUN git clone git://git.savannah.gnu.org/gsl.git && \
	cd /gsl && \
	./autogen.sh && \
	./configure && \
	make && \
	make install && \
	cd /
	
#Install htslib
RUN git clone https://github.com/samtools/htslib.git && cd htslib && make && cd /

#Install Samtools
RUN git clone git://github.com/samtools/samtools.git && cd samtools && autoheader && autoconf -Wno-syntax && ./configure && make && make install && cd /

#Install bcftools
RUN git clone git://github.com/samtools/bcftools.git && cd bcftools && make && cd /
#RUN git clone git://github.com/samtools/bcftools.git && cd bcftools && autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters && make && cd /

#Install Picard
RUN git clone https://github.com/broadinstitute/picard.git && cd picard/ && ./gradlew shadowJar && ./gradlew clean && cd /

#Install GATK4
RUN git clone https://github.com/broadinstitute/gatk.git && cd gatk/ && ./gradlew bundle && ./gradlew clean && cd /

#Install Java11
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 40976EAF437D05B5
#RUN apt install -y software-properties-common && \
#	apt-get update -q && \
#	add-apt-repository ppa:openjdk-r/ppa && \ 
#	apt-get update -q && \
#	apt install -y openjdk-11-jdk
#RUN apt install -y default-jdk
#RUN apt install -y software-properties-common && \
#	apt --allow-releaseinfo-change update && \
#	add-apt-repository ppa:webupd8team/java &&\
#	apt --allow-releaseinfo-change update && \
#	apt install -y dirmngr && \
#	apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EA8CACC073C3DB2A && \
#	echo "deb http://ppa.launchpad.net/linuxuprising/java/ubuntu bionic main" | tee /etc/apt/sources.list.d/linuxuprising-java.list && \
#	apt --allow-releaseinfo-change update && \
#	apt install -y oracle-java11-installer && \
#	apt install -y oracle-java11-set-default
#
#RUN echo "deb http://ftp.de.debian.org/debian stretch main" | tee /etc/apt/sources.list\ && \
#	apt update && \
#	apt install -y openjdk-11-jdk
#RUN apt install -y software-properties-common && \
#	add-apt-repository ppa:linuxuprising/java &&\
#	apt-get update && \
#	apt install -y --allow-unauthenticated oracle-java11-installer-local && \
#	apt install oracle-java11-set-default-local
#Move to start if it fixes IGV


#Install IGV Possible Debug. My test VM crashes at ./gradlew test but no errors thrown, so assumed working.
# RUN git clone https://github.com/igvteam/igv.git && cd igv/ && ./gradlew createDist && ./gradlew createToolsDist && ./gradlew test --no-daemon && cd /
#Skipped for now until java 11 figured out

#Install Trimmomatic
RUN git clone https://github.com/timflutre/trimmomatic.git && cd trimmomatic/ && make && make check && make install && cd /

#Install FASTQC
RUN git clone https://github.com/s-andrews/FastQC.git && cd FastQC/ && chmod 755 fastqc && ln -s /FastQC/fastqc /usr/local/bin/fastqc && cd /

#Install SPAdes
RUN curl -s "https://api.github.com/repos/ablab/spades/releases" | grep download | grep Linux.tar.gz | head -n 1 | awk '{print $2}' | xargs curl -L -o /SPAdes.tar.gz
RUN tar -zxvf SPAdes.tar.gz && rm SPAdes.tar.gz && mv SPAdes* SPAdes && cd /
ENV PATH $PATH:/SPAdes/bin/ 

#Install Freebayes
RUN git clone --recursive git://github.com/ekg/freebayes.git && \
	cd freebayes/ && \
	make && \
	make install && \
	cd /

#Install yaggo (prerequisite for Mummer)
RUN gem install yaggo

#Install Mummer
RUN curl -s "https://api.github.com/repos/mummer4/mummer/releases" | grep download | grep tar.gz | head -n 1 | awk '{print $2}' | xargs curl -L -o /mummer.tar.gz
RUN tar -zxvf mummer.tar.gz && rm mummer.tar.gz && mv mummer* mummer && \
	cd mummer/ && \
	./configure --prefix=/mummer/ && \
	make && \
	make install && \
	cd /
ENV PATH $PATH:/mummer/

#Install Prodigal (Prerequisite for Circlator)
RUN git clone https://github.com/hyattpd/Prodigal.git && cd Prodigal && make install && cd /

#Install Kraken2
RUN git clone https://github.com/DerrickWood/kraken2.git && \
	cd kraken2/ && \
	./install_kraken2.sh /kraken2/ && \
	cd /
ENV PATH $PATH:/kraken2/:/kraken2/kraken2-build/:/kraken2/kraken2-inspect/ \

#Install beast 1.x
RUN wget -O /beast1.tgz https://github.com/beast-dev/beast-mcmc/archive/v1.10.4.tar.gz && \
	tar -zxf beast1.tgz && \
	rm beast1.tgz && \
	mv beast-mcmc* beast1 && \
	cd /
ENV PATH $PATH:/beast1/bin/

#Install beast 2.x
#RUN curl -s "https://api.github.com/repos/CompEvol/beast2/releases/latest" | grep download | grep tgz | head -n 1 | awk '{print $2}' | xargs curl -L -o /beast2.tgz
#RUN tar -zxvf beast2.tgz && rm beast2.tgz && mv beast beast2 && export PATH=$PATH:/beast2/bin/ && cd /
#Temporary workaround, needs updating issue opened https://github.com/CompEvol/beast2/issues/875
RUN wget -O /beast2.tgz https://github.com/CompEvol/beast2/releases/download/v2.6.0/BEAST.v2.6.0.Linux.tgz && \
	tar -zxvf beast2.tgz && \
	rm beast2.tgz && \
	mv beast beast2 && \
	cd /
ENV PATH $PATH:/beast2/bin/ && \

#Install Figtree
RUN curl -s "https://api.github.com/repos/rambaut/figtree/releases/latest" | jq --arg PLATFORM_ARCH "tgz" -r '.assets[] | select(.name | endswith($PLATFORM_ARCH)).browser_download_url' | xargs curl -L -o /figtree.tgz
RUN tar -zxvf figtree.tgz && \
	rm figtree.tgz && \
	mv FigTree* FigTree && \
	cd /
ENV PATH $PATH:/FigTree/bin/ 

#Install SnpEff (SnpSift included)
RUN wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip && rm snpEff_latest_core.zip 
ENV PATH $PATH:/snpEff/
#Optional: The database directory can be changed in snpEff.config. Default is in the installation folder (./data/).
#See: http://snpeff.sourceforge.net for more information

#Install Circos #UI#
RUN wget -O /Circos.tgz http://circos.ca/distribution/circos-current.tgz && \
	tar -zxvf Circos.tgz && \
	rm Circos.tgz

#Install Miniconda (Prerequisite for MTBseq)
RUN wget -O /miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash /miniconda.sh -b -f -p /miniconda/ && rm /miniconda.sh && export PATH=$PATH:/miniconda/bin/ && \
	conda install -y anaconda && \
	conda install -y -c bioconda mtbseq && \
	mkdir /miniconda/mtbdependencies/ && \ 
	wget -O /miniconda/mtbdependencies/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 -U "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36" --referer https://software.broadinstitute.org/gatk/download/archive 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef' && \
	gatk3-register /miniconda/mtbdependencies/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
	conda install -y -c bioconda mykrobe && \
	conda install -y -c bioconda tb-profiler && \
	conda install -y -c bioconda pilon && \
	#conda install -y -c bioconda quast && \
	cd /

#Install QUAST alternate
RUN git clone https://github.com/ablab/quast.git && \
	cd quast && \
	./setup.py install_full &&\
	cd /

#Install VCFtools
RUN git clone https://github.com/vcftools/vcftools.git && \
	cd vcftools && \
	./autogen.sh && \
	./configure && \
	make && \
	make install && \
	cd /
	
#Install trimAl
RUN git clone https://github.com/scapella/trimal.git && \
	cd trimal/source && \
	make && \
	cd /
ENV PATH $PATH:/trimal/source/
#Install Racon
RUN git clone --recursive https://github.com/isovic/racon.git racon && \
	cd racon && \
	mkdir build && \
	cd build && \
	cmake -DCMAKE_BUILD_TYPE=Release .. && \
	make && \
	make install && \
	cd /

#Install Canu
RUN curl -s "https://api.github.com/repos/marbl/canu/releases/latest" | jq --arg PLATFORM_ARCH "Linux-amd64.tar.xz" -r '.assets[] | select(.name | endswith($PLATFORM_ARCH)).browser_download_url' | xargs curl -L -o /canu.tar.xz
RUN tar -xJf canu.tar.xz && \
	rm canu.tar.xz && \
	mv canu-* canu && \
	cd /
ENV PATH $PATH:/canu/Linux-amd64/bin

#Install Circlator
RUN pip3 install --upgrade setuptools
RUN pip3 install circlator && cd /

#Install TempEst
RUN wget -O /tempest.tgz 'http://tree.bio.ed.ac.uk/download.php?id=102&num=3'
RUN tar -zxvf tempest.tgz && \
	rm tempest.tgz && \
	mv TempEst* TempEst && \
	cd /
ENV PATH $PATH:/TempEst/bin/

#Install MEGA5 #May become outdated
RUN wget https://www.megasoftware.net/do_force_download/megax_10.0.5-1_amd64.deb #GUI version
RUN dpkg -i megax_10.0.5-1_amd64.deb && cd /

###GUI Installation
#Install and configure VNC server
RUN apt-get update && apt-get install -y tightvncserver twm
ENV DISPLAY :0
ENV USER root
RUN mkdir /root/.vnc
RUN echo password | vncpasswd -f > /root/.vnc/passwd
RUN chmod 600 /root/.vnc/passwd

# Install BEAGLE
# use latest release v3.0.1, issue #786
RUN git clone --branch v3.0.1 --depth=1 https://github.com/beagle-dev/beagle-lib.git
RUN cd /beagle-lib && ./autogen.sh && ./configure --prefix=/usr/local && make install
RUN ldconfig

# Ant build fails if the repo dir isn't named beast2
#RUN mkdir /root/beast2
#WORKDIR /beast2

#ADD . ./

CMD vncserver $DISPLAY -geometry 1920x1080; ant travis

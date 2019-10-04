#2-in-1 version, default.
#Based on Rocker's R-verse docker, with google cloud docker installed on top.
#To run the docker: docker run -it -p 59000:5900 qimrrebase /bin/bash
#From within the docker: vncserver $DISPLAY -geometry 1920x1080
#From a vnc client (e.g. Remmina), connect via VNC protocol to: localhost:5900
#To switch java versions: update-alternatives --config java

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


EXPOSE 80 443 22 9418 5900
#22 for SSH, 9418 for git, 80 and 443 general catch-alls, 5900 for VNC (GUI)

#Uncomment if deployment scenario does not have access to port 9418 which prevents git clones. This is a workaround.
#RUN git config --global url.https://github.com/.insteadOf git://github.com/

#Install apt packages, mostly dependencies, some are software
ENV IMAGE_PACKAGES="zlib1g-dev libz-dev libbz2-dev liblzma-dev libperl-dev libcurl4-gnutls-dev libgconf-2-4 libssl-dev libncurses5-dev libopenblas-base libtool libx11-dev libxext-dev libxrender-dev libxrandr-dev libxtst-dev libxt-dev libcups2-dev libasound2-dev pkg-config build-essential perl python-pip python3-pip autoconf automake jq ruby apache2 bwa gzip kalign tar wget vim ant bedtools"
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
RUN git clone git://github.com/samtools/samtools.git && cd samtools && \ 
	autoheader && autoconf -Wno-syntax && ./configure && make && make install && \
	cd /

#Install bcftools
RUN git clone git://github.com/samtools/bcftools.git && cd bcftools && make && cd /

#Install Picard
RUN git clone https://github.com/broadinstitute/picard.git && cd picard/ && ./gradlew shadowJar && ./gradlew clean && cd /

#Install GATK4
RUN git clone https://github.com/broadinstitute/gatk.git && cd gatk/ && ./gradlew bundle && ./gradlew clean && cd /

#Install Trimmomatic
RUN git clone https://github.com/timflutre/trimmomatic.git && cd trimmomatic/ && make && make check && make install && cd /

#Install FASTQC
RUN git clone https://github.com/s-andrews/FastQC.git && cd FastQC/ && chmod 755 fastqc && \
	ln -s /FastQC/fastqc /usr/local/bin/fastqc && \
	cd /

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
ENV PATH $PATH:/kraken2/:/kraken2/kraken2-build/:/kraken2/kraken2-inspect/

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
ENV PATH $PATH:/beast2/bin/

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
RUN wget -O /miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	bash /miniconda.sh -b -f -p /miniconda/ && \
	rm /miniconda.sh
ENV PATH $PATH:/miniconda/bin/
RUN conda install -y anaconda && \
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

#Install MEGA5 #May become outdated MOVE libgtk2 TO FRONT IF THIS WORKS
RUN wget https://www.megasoftware.net/do_force_download/megax_10.0.5-1_amd64.deb #GUI version
RUN apt install -y libgtk2.0-bin
RUN dpkg -i megax_10.0.5-1_amd64.deb && cd /

#Install Java11 (downloads prebuilt Java10 first to enable building of Java11)
RUN curl -s 'https://api.adoptopenjdk.net/v2/info/releases/openjdk10' | grep binary_link | grep OpenJDK10_x64_Linux | head -n 1 |  awk '{print $2}' | xargs echo | rev | cut -c 2- | rev | xargs curl -L -o /OpenJDK10u.tar.gz
RUN tar -zxvf OpenJDK10u.tar.gz && \
	mv jdk-10* jdk10 && \
	rm OpenJDK10u.tar.gz && \
	chmod a+x /jdk10/bin/java
ENV JAVA_HOME=/jdk10 BOOT_JDK=/jdk10 PATH=$JAVA_HOME/bin:$PATH
#J11 builder dependencies, move to start if works.
RUN apt install -y 
RUN git clone https://github.com/AdoptOpenJDK/openjdk-build.git && \
	cd /openjdk-build && \
	./makejdk-any-platform.sh -J /jdk10 --configure-args --disable-warnings-as-errors jdk11u && \
	cd /
	
#Install IGV Possible Debug. My test VM crashes at ./gradlew test but no errors thrown, so assumed working.
RUN git clone https://github.com/igvteam/igv.git && cd igv/ && ./gradlew createDist && ./gradlew test --no-daemon && cd /

###GUI Installation
#Install desktop environment
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt install -y xfce4 task-xfce-desktop

#Install and configure VNC server
RUN apt install -y tightvncserver twm
ENV DISPLAY=:0 USER=root
RUN mkdir /root/.vnc && \
	echo password | vncpasswd -f > /root/.vnc/passwd
COPY xstartup /root/.vnc/
RUN chmod 600 /root/.vnc/passwd && \
	chmod a+x /root/.vnc/xstartup

# qimr_mtb
Pipeline for mtb analysis built for QIMR. 
Please report bugs in Issues, some software are undergoing transitions from Java8 to Java11, which requires updating of the Dockerfile. 

## Manual Instructions
### Written for debian-based OS. May work on similar architectures with no or minor syntax tweaks.
#### Pull cloud build. (Not always reliable, cloud build tends to fail for such a large image.=)
1. Set up and log in to your docker account if you haven't https://www.docker.com/
2. Pull the image with `docker pull dockersubtest/qimr_mtb`
3. Make a shared directory, such as `mkdir ~/dhost_mount`
4. Run the docker with `docker run -it -p 5900:5900 -p 8787:8787 -v $HOME/dhost_mount:/dcont_mount/ qimr_mtb /bin/bash`

#### Local build
1. Set up and log in to your docker account if you haven't https://www.docker.com/
2. Ensure `git` is installed
3. Clone this git project `git clone https://github.com/Firedrops/qimr_mtb.git` and navigate to the folder (usually ~/qimr_mtb/)
4. Build the image `sudo docker build -t qimr_mtb .`
5. Wait... build is expected to take about 2 hours with a docker image ~40 GB.
6. Make a shared directory, such as `mkdir ~/dhost_mount`
7. Run the docker with `docker run -it -p 5900:5900 -p 8787:8787 -v $HOME/dhost_mount:/dcont_mount/ qimr_mtb /bin/bash`

### To use the GUI
1. Ensure you have a VNC client installed on the host machine (e.g. Remmina)
2. Start the VNC server from within the docker `vncserver $DISPLAY -geometry 1920x1080`
3. In your VNC client, select VNC protocol and connect to `localhost:5900` The password is simply `password`.

### To switch Java versions
In the docker, use `update-alternatives --config java`

### To exit the container (e.g. shut it down)
Simply enter `exit` command or `ctrl + c`.

### To leave the container without stopping it
#### Analogous to tabbing out to another window and letting it continue running in the background.
1. `ctrl + p` or `ctrl + q`
2. To attach back, `docker ps` (may have to add `sudo`), note the first 3 characters under `ID`. For example, JHD898A. 2 will be sufficient if you only have 1 or a few docker containers running, without ID overlaps.
3. `docker attach <ID>`. For example `docker attach JHD`.
#### This is not recommended. Consider simply by running another terminal window on the host machine if possible.

### To pass files in/out from host to docker
This can be useful for easily passing scripts, source files, and output files in/out of the docker.
Simply move/copy them into `~/dhost_mount` (host-side), and `/dcont_mount/` (docker-side). Note that the host is in home, and docker in root.
#### IMPORTANT! Only files stored to `/dcont_mount/` directory in the docker will be persistent. All other files will be lost on termination of docker.

### To use R (deprecated, can be used via VNC, but this method should still work)
1. Run any internet browser on your host machine.
2. Navigate to `http://localhost:8787`.
Further documentation can be found [here](https://ropenscilabs.github.io/r-docker-tutorial/02-Launching-Docker.html)

## Auto Instructions (Not yet functional)
### Written for debian-based OS. May work on similar architectures with no or minor syntax tweaks.
1. Set up a docker account if you haven't https://www.docker.com/
2. Ensure `git` is installed
3. Clone this git project `git clone https://github.com/Firedrops/qimr_mtb.git` and navigate to the folder (usually ~/qimr_mtb/)
4. Open `duser.txt`, `dpass.txt`, and `lpass.txt`. Replace `placeholder` text with your docker username, docker password, and linux password, respectively.
5. Optional: In case your working directory is not `~/qimr_mtb/`, open and edit the first few lines of `init.sh` to match your working directory.
6. Run `init.sh`.
  init.sh checks for existing docker installation, and installs it automatically if it is not found, and
    logs in to docker, and
      pulls the qimr_mtb docker image, and
        creates a directory `~/dhost_mount` to share files in/out of the docker, and
          runs and and attaches user input into the docker image.
  At this point (you can see when your user input becomes `root@xx...xx`) you are inside the docker, and can use its tools.

## Tools included in this docker
Bcftools

Beast 1 and 2

Bedtools

Bwamem 

Canu 

Circos

Circulator

FastQC

Figtree

Freebayes

Gatk

IGV

Kraken

MEGA

MTBseq

Minimap2

Mummer (and Yaggo)

Mykrobe predictor

## Citation
If you have found this useful in your research, please consider citing the authors: Larry Cai, Arnold Bainomugisa, and Lachlan Coin.
Picardtools

Pilon

Quast

R/R-Studio

RAxML

Racon

Rapiddr

Samtools

SnpEFF/snpSIFT

Spades

TempEst

Trimal

Trimmomatic

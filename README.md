# qimr_mtb
Pipeline for mtb analysis built for QIMR. 

## Instructions
### Written for debian-based OS. May work on similar architectures with small or minor syntax tweaks.
1. Set up a docker account if you haven't https://www.docker.com/
2. Ensure up to date repository and `git` is installed
3. `git clone https://github.com/Firedrops/qimr_mtb.git` and navigate to the folder (usually ~/qimr_mtb/)
4. Open `duser.txt`, `dpass.txt`, and `lpass.txt`. Replace `placeholder` text with your docker username, docker password, and linux password, respectively.
4.1 Optional: In case your working directory is not `~/qimr_mtb/`, open and edit the first few lines of `init.sh` to match your working directory.
5. Run `init.sh`.
  init.sh checks for existing docker installation, and installs it automatically if it is not found, and
    logs in to docker, and
      pulls the qimr_mtb docker image, and
        creates a directory `~/dhost_mount` to share files in/out of the docker, and
          runs and and attaches user input into the docker image.
  At this point (you can see when your user input becomes `root@xx...xx`) you are inside the docker, and can use its tools.

## To exit the container (e.g. shut it down)
Simply enter `exit` command or `ctrl + c`.

## To leave the container without stopping it
### This is analogous to tabbing out to another window and letting it continue running in the background
### Note: This is not strictly necessary. You can interface with the host machine simply by running another terminal window.
1. `ctrl + p` or `ctrl + q`
2. To attach back, `docker ps` (may have to add `sudo`), note the first 3 characters under `ID`. For example, JHD898A. 2 will be sufficient if you only have 1 or a few docker containers running, without ID overlaps.
3. `docker attach <ID>`. For example `docker attach JHD`.

## To pass files in/out from host to docker
This can be useful for easily passing scripts, source files, and output files in/out of the docker.
Simply move/copy them into `~/dhost_mount` (host-side), and `/dcont_mount/` (docker-side). Note that the host is in home, and docker in root.
## IMPORTANT! Only files stored to `/dcont_mount/` directory in the docker will be persistent. All other files will be lost on termination of docker.

## To use R
1. Run any internet browser on your host machine.
2. Navigate to `http://localhost:8787`.
Further documentation can be found [here](https://ropenscilabs.github.io/r-docker-tutorial/02-Launching-Docker.html)

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

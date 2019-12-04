1. Pipe 1 example usage:```docker run -it --rm -v /Backup/Data/Projects/PNG_TB/:/data/ dockersubtest/qimr_mtb:slim /data/ RB17MT1804 /data/H37Rv_refe /data/temp /data/out```
  Replace file name, reference path (exclude ".fasta"), and designated temp/out folders as required.
  Alternatively, on slurm:
  ```
sbatch --array 1-n /home/larry/qimr_mtb/Slim/runpipe.sh
  ```
2. Pipe 2 example usage: ```docker run -it --rm --entrypoint /pipe2.sh -v /home/lcoin/nextflow/data:/data/ -v /home/lcoin/nextflow/out:/out/ qimr_slim /data/ /data/H37Rv_refe /out/temp /out/out```
  Alternatively, on slurm:
  ```
sbatch --array 1-n /home/larry/qimr_mtb/Slim/runpipe2.sh
  ```
3. Using the above example, for each `.fasta`, there will be 2 unique folders generated, `./temp_RB17MT1804` and `./out_RB17MT1804`, automatically named after the input file name to avoid conflicts. These will then be used in `pipe2.sh`

#Refreshing Docker on slurm:
Due to some bugs with dockerhub, pulling a build from dockerhub results in errors.
1. Commit on github, trigger manual build on dockerhub.
2. Build on own PC/nectar `docker build -t qimr_slim .`
3. Run an interactive container in detached mode `docker run -itd --entrypoint /bin/bash qimr_slim`
4. Check for container id `docker container ls`
5. Commit built image to dockerhub image `docker commit ae83b dockersubtest/qimr_mtb:slim`
6. Push built image `docker push dockersubtest/qimr_mtb:slim`
7. Spin down docker container `docker stop ae8eb`

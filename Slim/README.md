# Usage Instructions
1. Pipe 1 example usage:```docker run -it --rm -v /Backup/Data/Projects/PNG_TB/:/data/ dockersubtest/qmrl_mtb:slim /data/ RB17MT1804 /data/H37Rv_refe /data/temp /data/out```
  Replace file name, reference path (exclude ".fasta"), and designated temp/out folders as required.
  Alternatively, on slurm:
  ```
sbatch --array 1-n /home/larry/qmrl_mtb/Slim/runpipe.sh
  ```
2. Pipe 2 example usage: ```docker run -it --rm --entrypoint /pipe2.sh -v /home/lcoin/nextflow/data:/data/ -v /home/lcoin/nextflow/out:/out/ qmrl_slim /data/ /data/H37Rv_refe /out/temp /out/out```
  Alternatively, on slurm:
  ```
sbatch --array 1-n /home/larry/qmrl_mtb/Slim/runpipe2.sh
  ```
3. Using the above example, for each `.fasta`, there will be 2 unique folders generated, `./temp_RB17MT1804` and `./out_RB17MT1804`, automatically named after the input file name to avoid conflicts. These will then be used in `pipe2.sh`

# Refreshing Docker on slurm:
Due to some bugs with dockerhub, pulling a build from dockerhub results in errors.

Auto: run `fullupdate.sh` in root git folder.

Manual:
1. Commit on github, trigger manual build on dockerhub.
2. Navigate to `<home>/qmrl_mtb/Slim` folder OR adjust below path arguments accordingly.
3. Build on own PC/nectar `docker build -t qmrl_slim .`
4. Run an interactive container in detached mode `docker run -itd --entrypoint /bin/bash qmrl_slim`
5. Check for container ID `docker container ls`. For example, it is `ae83ba82s`.
6. Commit built image to dockerhub image `docker commit ae83 dockersubtest/qmrl_mtb:slim`
7. Push built image `docker push dockersubtest/qmrl_mtb:slim`
8. Spin down docker container `docker stop ae8eb`
9. Clear legacy containers on slurm `sbatch --array 1-5 clearslurm.sh`

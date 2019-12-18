#!/bin/bash
git pull
cd Slim/
docker build -t qimr_slim .
docker run -itd --entrypoint /bin/bash qimr_slim
curr_docker=`docker container ls | grep qimr | awk '{print $1}'`
docker commit $curr_docker dockersubtest/qimr_mtb:slim
docker push dockersubtest/qimr_mtb:slim
docker stop $curr_docker
unset curr_docker
sbatch --array 1-5 clearslurm.sh
cd ..

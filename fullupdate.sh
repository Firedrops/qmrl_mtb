#!/bin/bash
git pull
cd Slim/
docker build -t qmrl_slim .
docker run -itd --entrypoint /bin/bash qmrl_slim
curr_docker=`docker container ls | grep qmrl | awk '{print $1}'`
docker commit $curr_docker dockersubtest/qmrl_mtb:slim
docker push dockersubtest/qmrl_mtb:slim
docker stop $curr_docker
unset curr_docker
sbatch --array 1-5 clearslurm.sh
cd ..

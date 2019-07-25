export USB="insert usb path here"

#portable docker installer



#save pre-built docker image for offline portability
sudo docker save -o $USB/docker/qimr_mtb.docker qimr_mtb

#To load docker image on the new machine, run this command
#sudo docker load -i qimr_mtb.docker
#make sure the working directory is in the where the .docker file is, or just add the full directory,
#for example, /home/<username>/qimr_mtb.docker

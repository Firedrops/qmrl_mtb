export USB="insert usb path here"
#retain the double quotes

#save pre-built docker image for offline portability
sudo docker save -o $USB/docker/qmrl_mtb.docker qmrl_mtb

#To load docker image on the new machine, run this command
#sudo docker load -i qmrl_mtb.docker
#make sure the working directory is in the where the .docker file is, or just add the full directory,
#for example, /home/<username>/qmrl_mtb.docker

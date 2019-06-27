#!/bin/bash

export DockerUsername="foo" #replace foo in quotes with your username, retain the quotation marks.
export lpass="~/qimr_mtb/lpass.txt" #path to .txt containing linux password
#By default, sudo password is only prompted once every 15 minutes, so I have chosen to pass it only once per if loop to simplify commands.
export dpass="~/qimr_mtb/dpass.txt" #path to .txt containing docker password

#Checks if Docker is installed
echo "Checking if Docker is installed"
dockercheck1=$(docker run hello world | grep "working correctly" | wc -l)
if [ "$dockercheck1 -ge 1" ]; then
  echo "Docker appears to be installed, proceeding."
else
  echo "Docker not found, attempting installation."
  cat $lpass | sudo -S apt-get -y update
  sudo apt-get -y install apt-transport-https ca-certificates curl gnupg-agent software-properties-common
  curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -
  sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
  sudo apt-get -y update
  sudo apt-get -y install docker-ce docker-ce-cli containerd.io
  echo "sudo apt-get -y install apt-transport-https ca-certificates curl gnupg-agent software-properties-common"
  echo "curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -"
  echo "sudo apt-get -y update"
  echo "sudo apt-get -y install docker-ce docker-ce-cli containerd.io"
  #Checks for successful installation
  dockercheck2=$(sudo docker run hello-world | grep "working correctly" | wc -l)
  if [ "$dockercheck2 -ge 1" ]; then
    echo "Docker has been successfully installed, proceeding."
  else
      echo "Docker could not be installed. Aborting."
      exit 1
  fi
fi

#Log in to docker using username declared at the top and password from dpass.txt
echo "Logging in to Docker with supplied credentials."
dockerlogin=$(cat $dpass | docker login --username $DockerUsername --password-stdin | grep "Succeeded" | wc -l)
if [ "$dockerlogin -ge 1"]; then
  echo "Logged in to Docker successfully."
else
  echo "Log in failed. Aborting. Please doublecheck declared username and dpass.txt. Case-sensitive."
  exit 1
fi

#Pull Docker image #WIP
docker pull dockersubtest/qimr_mtb:latest
docker run -itd qimr_mtb

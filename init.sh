export DockerUsername="foo" #replace foo in quotes with your username, retain the quotation marks.

#Checks if Docker is installed
dockercheck1=$(docker run hello world | grep "working correctly" | wc -l)
if [ "$dockercheck1 -ge 1" ]; then
  echo "Docker appears to be installed, proceeding."
else
  echo "Docker not found, attempting automatic installation."
  sudo apt-get -y update
  sudo apt-get -y install docker-ce docker-ce-cli containerd.io
  echo "sudo apt-get -y update"
  echo "sudo apt-get -y install docker-ce docker-ce-cli containerd.io"
  #Checks for successful installation
  dockercheck2=$(docker run hello world | grep 'working\|correctly' | wc -l)
  if [ "$dockercheck2 -ge 1" ]; then
    echo "Docker has been successfully installed, proceeding."
  else
      echo "Docker could not be installed. Aborting."
      exit 1
  fi
fi

#Logging in to docker using username declared at the top and password from my_password.txt
dockerlogin=$(cat ~/my_password.txt | docker login --username $DockerUsername --password-stdin | grep "Succeeded" | wc -l)
if [ "$dockerlogin -ge 1"]; then
  echo "Logged in to Docker successfully."
else
  echo "Log in failed. Aborting. Please doublecheck declared username and my_password.txt. Case-sensitive."
  exit 1
fi

#Pull Docker image #WIP
docker pull dockersubtest/qimr_mtb:latest
docker run -itd qimr_mtb

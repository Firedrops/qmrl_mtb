#!/bin/bash

#SBATCH --job-name=beast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16000 # mb
#SBATCH --time=9999:00:00

JAVA_LIB="/swold/beast/BEASTv1.8.2/lib:/usr/local/lib:/home/lcoin/lib"
JAVA_CP="/swold/beast/BEASTv1.8.2/lib/beast.jar:/swold/beast/BEASTv1.8.2/lib/beast-beagle.jar"
#JAVA_LIB="/swold/beast/BEASTv1.10.2/lib:/usr/local/lib:/home/lcoin/lib:"
#JAVA_CP="/swold/beast/BEASTv1.10.2/lib/beast.jar:/swold/beast/BEASTv1.10.2/lib/beast-beagle.jar"

#java -Xms64m -Xmx8092m -Djava.library.path=/swold/beast/BEASTv1.10.2/lib:/usr/local/lib:/home/lcoin/lib: -cp /swold/beast/BEASTv2.10.2/lib/beast.jar:/swold/beast/BEASTv1.10.2/lib/beast-beagle.jar dr.app.beast.BeastMain DARU_95_real_f5_beast1_trial1_3.xml
java -Xms64m -Xmx8092m -Djava.library.path=$JAVA_LIB -cp $JAVA_CP dr.app.beast.BeastMain DARU_95_real_f5_beast1_trial1_3.xml 

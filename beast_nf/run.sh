#!/bin/bash

#SBATCH --job-name=Kallisto
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64000 # mb
#SBATCH --time=9999:00:00

JAVA_LIB="/swold/beast/BEASTv1.8.2/lib:/usr/local/lib:/swold/beast/beagle-lib/lib/"
JAVA_CP="/swold/beast/BEASTv1.8.2/lib/beast.jar:/swold/beast/BEASTv1.8.2/lib/beast-beagle.jar"

java -Xms64m -Xmx8092m -Djava.library.path=$JAVA_LIB -cp $JAVA_CP dr.app.beast.BeastMain ./data/DARU_95_real_f5_beast1_trial1_3.xml 

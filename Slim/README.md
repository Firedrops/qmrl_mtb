1. Install and log in to docker
2. Example usage: ```docker run -it --rm -v /Backup/Data/Projects/PNG_TB/:/data/ qimr_slim /data/ RB17MT1804 /data/H37Rv_refe /data/temp /data/out```
  Replace file name, reference path (exclude ".fasta"), and designated temp/out folders as required.
  Alternatively, 
  ```
ls /data/ | grep _R1.fastq > filesin.txt

n=$(wc -l filesin.txt)
sbatch --array 1-$n runpipe.sh

  ```

3. Using the above example, there will be 2 unique folders generated, "./temp_RB17MT1804" and "./out_RB17MT1804", automatically named after the input file name to avoid conflicts.

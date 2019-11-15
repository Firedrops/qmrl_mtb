1. Install and log in to docker
2. Example usage: ```docker run -it --rm -v /Backup/Data/Projects/PNG_TB/:/data/ qimr_slim /data/ RB17MT1804 /data/H37Rv_refe /data/temp /data/out```
  Replace file name, reference path (exclude ".fasta"), and designated temp/out folders as required.
  Alternatively, 
  ```
ls | grep _R1.fastq | cut -d _ -f 1 > /home/larry/png_mtb/filesin.txt
n=$(wc -l /home/larry/png_mtb/filesin.txt)
sbatch --array 1-3 /home/larry/qimr_mtb/Slim/runpipe.sh

  ```

3. Using the above example, there will be 2 unique folders generated, "./temp_RB17MT1804" and "./out_RB17MT1804", automatically named after the input file name to avoid conflicts.

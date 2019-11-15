1. Install and log in to docker
2. Example usage: ```docker run -it --rm -v /Backup/Data/Projects/PNG_TB/:/data/ qimr_slim /data/ RB17MT1804 /data/H37Rv_refe /data/temp /data/out```
  Replace file name, reference path (exclude ".fasta"), and designated temp/out folders as required.
  Alternatively, 
  ```
  for i in /data/; 
    do docker run -it --rm -v /Backup/Data/Projects/PNG_TB/:/data/ qimr_slim /data/ $i /data/H37Rv_refe /data/temp /data/out
  done
  ```

3. Using the above example, there will be 2 unique folders generated, "./temp_RB17MT1804" and "./out_RB17MT1804", automatically named after the input file name to avoid conflicts.

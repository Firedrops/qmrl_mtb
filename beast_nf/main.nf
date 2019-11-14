#!/usr/bin/env nextflow



/*
 * Define the default parameters
 */ 

params.xml     = "$baseDir/data/DARU_95_real_f5_beast1_trial1_3.xml"

JAVA_LIB="/swold/beast/BEASTv1.8.2/lib:/usr/local/lib:/swold/beast/beagle-lib/lib/"
JAVA_CP="/swold/beast/BEASTv1.8.2/lib/beast.jar:/swold/beast/BEASTv1.8.2/lib/beast-beagle.jar"

xml_file     = file(params.xml)



/**********
 * Process 1: run beast
 *
 **********/
process 'run_beat' {
  tag "$xml.baseName"

  
  publishDir '/data/output/'
  
  input:
      file xml from xml_file
  output:
      file "${xml.baseName}.trees" into trees_ch
	  file "${xml.baseName}.log" into log_ch
	   file "${xml.baseName}.tress.txt" into trees_ch1

  script:
  """
  java -Xms64m -Xmx8092m -Djava.library.path=$JAVA_LIB -cp $JAVA_CP dr.app.beast.BeastMain 
  """
}
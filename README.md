# Microbiome-Analysis

Analysis of mibrobial DNA in dairy heifer pilot study. 

Final files to be passed into CLC Genomic Workbench. 

- `Stats_QC.py`: This script contains two functions which
  - Gets statistics on the fastq file passed into it and prints the number of reads, number of bases, average sequence length, and GC content. 
  - Finds and exports the quality and base content information on passed in files
  
- `Tables_Plots.R`: This script is designed to read in files created by `Stats_QC.py` and 
  - Create a table of statistics
  - Create a plot of quality scores
  - Create a plot of percentages of bases at specified locations on the reads  

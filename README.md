# Nanopore Pipeline
Pipeline for analysis and SV calling of targeted nanopore sequencing data.

## Usage 
nanopore_pipeline.sh -n \<run name\> -d \<directory containing fast5 files\> -o \<output directory\> 
  
  optional arguments for adaptive sampling:
  
  -b:  _bed file containing targeted regions_
  
  -a:  _adaptive sampling summary output file_
  

## Dependencies

Guppy,
minimap2,
Python 3+,
R + tidyverse,
pycoQC,
samtools,
bedtools,
mosdepth,
cuteSV.

# Nanopore Pipeline
Pipeline for basecalling, alignment and SV calling of targeted nanopore sequencing data from an adaptive sampling sequencing run.

## Usage 
adaptive_panel.sh -n *run name* -d *\path\to\sequencing\data* -b *\path\to\adaptive sampling.bed* -o *\path\to\output\directory*
                      -m *\path\to\reference.mmi* -r *\path\to\reference\fasta.fa* 

Optional: -a Option to skip adaptive sampling analysis Y/N. Default N


## Dependencies

Dorado with minimap2,
Sniffles2,
Python 3+,
R + tidyverse,
NanoPlot,
samtools,
bedtools,
mosdepth,
SvAnna,
cuteSV.

#!/bin/#!/usr/bin/env bash
#analysis of adaptive sampling adaptive_stats
#TO DO: change from bedtools to mosdepth

#source /home/nanopore/miniconda3/etc/profile.d/conda.sh

#parse arguments
while getopts d:n:s:b:w: flag; do
  case "${flag}" in
    d) pipeline_dir="${OPTARG}";;
    n) run_name="${OPTARG}";;
    s) adaptive_summary="${OPTARG}";;
    b) bed_file="${OPTARG}";;
    w) work_dir="${OPTARG}";;
  esac
done


mkdir -p "$work_dir"/adaptive_stats/adaptive_coverage
cd "$work_dir"/adaptive_stats


#create bam files containg read ids for each adaptive sampling decision
echo "Subseting bam file..."
# -b bam file
# -a adaptive sequencing summary file
#
python3 $pipeline_dir/SCRIPTS/extract_reads_adaptive.py \
--bam "$work_dir"/alignment/"$run_name".bam \
--adaptive_output "$adaptive_summary" \
--out "$run_name".bam

#get list of bam files from last step
ls *.bam > bam_files.txt

#index subsetted files
while read bam_file; do
  samtools sort -@ 20 -o $bam_file.sorted.bam $bam_file
  samtools index $bam_file.sorted.bam
done < bam_files.txt

#calculate on target percentages (bedtools)
#mkdir ./ADAPTIVE_COVERAGE

for f in *.sorted.bam; do
  #create directories for each adaptive decision
  decision=${f%%.bam*}
  mkdir ./adaptive_coverage/"$decision"
  #run bedtools to calculate coverage summary and per base depth of each feature in bed_file
  bedtools coverage -a $bed_file -b $f -d > ./adaptive_coverage/"$decision"/"$decision"_per_base_depth.tsv
done


#stats
#TO ADD TO R SCRIPT
##per gene coverage
##depth calculations


#descriptive stats from adaptive sampling
Rscript $pipeline_dir/SCRIPTS/adaptive_stats.r "$adaptive_summary" $run_name

cd ./ADAPTIVE_COVERAGE

#depth and coverage calculations on .tsv output from bedtools
#T
Rscript $pipeline_dir/SCRIPTS/coverage_adaptive_panel.r *.tsv $run_name

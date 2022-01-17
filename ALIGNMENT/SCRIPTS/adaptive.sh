#!/bin/#!/usr/bin/env bash
#analysis of adaptive sampling adaptive_stats

source /home/nanopore/miniconda3/etc/profile.d/conda.sh

#parse arguments
while getopts d:n:s:b: flag; do
  case "${flag}" in
    d) pipeline_dir="${OPTARG}";;
    n) run_name="${OPTARG}";;
    s) adaptive_summary="${OPTARG}";;
    b) bed_file="${OPTARG}";;
  esac
done


mkdir ./adaptive_stats
cd ./adaptive_stats

#create bam files containg read ids for each adaptive sampling decision
echo "Subseting bam file..."
python $pipeline_dir/SCRIPTS/extract_reads_adaptive.py -b ../alignment/"$run_name".bam -a $adaptive_summary -o "$run_name".bam

#get list of bam files from last step
ls *.bam > bam_files.txt

#index subsetted files
while read bam_file; do
  samtools sort -o $bam_file.sorted.bam $bam_file
  samtools index $bam_file.sorted.bam
done < bam_files.txt

#calculate on target percentages (bedtools)
#IN PROGRESS
mkdir ./COVERAGE

for f in *.sorted.bam; do
  #create directories for each adaptive decision
  decision=${f%%.bam*}
  mkdir ./COVERAGE/"$decision"
  #run bedtools to calculate depth of each feature in bed_file
  bedtools coverage -a $bed_file -b $f -d > ./COVERAGE/"$decision"/"$f"_per_base_depth.tsv
done

#stats
conda deactivate
conda activate r_env
Rscript $pipeline_dir/SCRIPTS/adaptive_stats.r $adaptive_summary $run_name

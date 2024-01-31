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

# echo "PIPELINE DIR: $pipeline_dir"
# echo "RUN NAME: $run_name"
# echo "ADAPTIVE SUMMARY: $adaptive_summary"
# echo "BED FILE: $bed_file"
# echo "WORK DIR: $work_dir"
# echo "BAM FILE: "$work_dir"/alignment/"$run_name".bam"

mkdir -p "$work_dir"/adaptive_stats/depth
#cd "$work_dir"/adaptive_stats


#create bam files containg read ids for each adaptive sampling decision
echo $(date)
echo "INFO: Subseting bam file..."

###### subset 
python3 $pipeline_dir/SCRIPTS/subset_adaptive.py \
  --adaptive_output "$adaptive_summary" \
  --output_dir "$work_dir"/adaptive_stats \
  --run_name "$run_name"

#subset bam file using samtools
samtools view \
  -@ 16 \
  -hN "$work_dir"/adaptive_stats/"$run_name"_stop_receiving_read_ids.txt \
  "$work_dir"/alignment/"$run_name".bam > "$work_dir"/alignment/"$run_name"_stop_receiving.bam

samtools sort \
  -@16 \
  -o "$work_dir"/alignment/"$run_name"_stop_receiving.sorted.bam \
  "$work_dir"/alignment/"$run_name"_stop_receiving.bam

samtools index \
  "$work_dir"/alignment/"$run_name"_stop_receiving.sorted.bam

#TODO bedtools coverage for stop receiving bam file
bedtools coverage \
  -a "$bed_file" \
  -b "$work_dir"/alignment/"$run_name"_stop_receiving.sorted.bam \
  -d > "$work_dir"/adaptive_stats/depth/"$run_name"_stop_receiving_per_base_depth.tsv

echo $(date)
echo "INFO: Running descriptive statistics..."
#descriptive stats from adaptive sampling
Rscript $pipeline_dir/SCRIPTS/adaptive_stats.r \
  "$adaptive_summary" \
  $run_name \
  "$work_dir"/adaptive_stats/depth

#cd ./ADAPTIVE_COVERAGE

#depth and coverage calculations on .tsv output from bedtools

Rscript $pipeline_dir/SCRIPTS/coverage_adaptive_panel.r \
  "$work_dir"/adaptive_stats/depth/"$run_name"_stop_receiving_per_base_depth.tsv\
  "$run_name" \
  "$work_dir"/adaptive_stats/depth

####

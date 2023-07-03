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

mkdir -p "$work_dir"/adaptive_stats/adaptive_coverage
#cd "$work_dir"/adaptive_stats


#create bam files containg read ids for each adaptive sampling decision
echo $(date)
echo "INFO: Subseting bam file..."
# -b bam file
# -a adaptive sequencing summary file
######
#TODO change this section to use samtools view -N on output of subset_adaptive.py
python3 $pipeline_dir/SCRIPTS/subset_adaptive.py \
--adaptive_output "$adaptive_summary" \
--output_dir "$work_dir"/adaptive_stats \
--run_name "$run_name"

#subset bam file using samtools
samtools view -@ 16 -hN "$work_dir"/adaptive_stats/"$run_name"_stop_receiving_read_ids.txt \
"$work_dir"/alignment/"$run_name".bam > "$work_dir"/alignment/"$run_name"_stop_receiving.bam

samtools sort -@16 -o "$work_dir"/alignment/"$run_name"_stop_receiving.sorted.bam "$work_dir"/alignment/"$run_name"_stop_receiving.bam
samtools index "$work_dir"/alignment/"$run_name"_stop_receiving.sorted.bam
exit 1

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

####

#stats
#TO ADD TO R SCRIPT
##per gene coverage
##depth calculations

echo $(date)
echo "INFO: Running descriptive statistics..."
#descriptive stats from adaptive sampling
Rscript $pipeline_dir/SCRIPTS/adaptive_stats.r "$adaptive_summary" $run_name

#cd ./ADAPTIVE_COVERAGE

#depth and coverage calculations on .tsv output from bedtools

Rscript $pipeline_dir/SCRIPTS/coverage_adaptive_panel.r ./adaptive_coverage/stop_receiving_"$run_name"/stop_receiving_"$run_name".tsv "$run_name" ./

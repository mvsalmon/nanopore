#!/bin/bash

#####
##
## Pipeline for analysis of targeted Nanopore sequencing data using adaptive sampling 
## Inteded for the detection of gene fusions in myeloid neoplasms
## Wessex genomics laboratory
##
####

#Usage
helpFunction()
{
   echo ""
   echo "Pipeline for analysis of Nanopore runs."
   echo ""
   echo "Usage: $0 -h display help -n run_name -d run_dir -o output_dir -b bed_file"
   echo "Required arguments:"
   echo -e "\t-n Name of Nanopore run"
   echo -e "\t-d Path to directory containing Nanopore run data. Should be the experiment level directory if several runs need to be analysed together."
   echo -e "\t-b BED file for adaptive sampling analysis" 
   echo -e "\t-o Output directory path"
   echo -e "\t-m Path to reference genome mmi index"
   echo -e "\t-r Path to reference genome fa file"
   echo -e "\t-a Adaptive sampling analysis is enabled by default. Set this flag to 0 to skip this analysis"
   echo -e "\t-p Phenotype term to pass to SvAnna for SV prioritisation."

   exit 1 # Exit script after printing help
}

#set defaults
adaptive_sampling=1

#parse arguments
#TODO change this from getopts
while getopts n:d:b:o:m:r:p:s:l:q:v:ha opt; do
  case "$opt" in
    n) run_name="$OPTARG";;
    d) run_dir="$OPTARG";;
    b) bed_file="$OPTARG";;
    o) output_dir="$OPTARG";;
    m) mmi_index="$OPTARG";;
    r) ref_index="$OPTARG";;
    p) phenotype_term="$OPTARG";;
    #options to skip steps - useful for isolating/debugging specific aspects
    s) skip_basecalling="$OPTARG";;
    l) skip_alignment="$OPTARG";;
    q) skip_qc="$OPTARG";;
    v) skip_SV="$OPTARG";;
    #p) skip_adaptive="$OPTARG";;
    h) helpFunction;;
    a) adaptive_sampling=0;;
    
  esac
done

#Help if mandatory arguments are empty
if [ -z "$run_name" ] || [ -z "$run_dir" ] || [ -z "$output_dir" ] || \
   [ -z "$mmi_index" ] || [ -z "$bed_file" ] || [ -z "$ref_index" ] || \
   [ -z "$phenotype_term" ]
   
then
   echo ""
   echo "ERROR: Missing one or more required arguments. See help message below.";
   helpFunction 
   exit 1
fi

# Check for running dorado service before use and exit if running 
# This frees up pre-allocated dorado resources for the analysis

pid=$( nvidia-smi | grep dorado | awk '{print $5}' )

if [[ "$pid" =~ ^[0-9]+$ ]]; then
  >&2 echo "EXITING: Running Dorado instance detected. Try: 'sudo service doradod stop' then retry."
  exit 1

 else
   >&2 echo $(date)
   >&2 echo "INFO: No running Dorado detected, running new analysis..."
fi

#####MAIN PIPELINE######

#create analysis dirs
#dir variables
pipeline_dir=$(pwd)
work_dir="$output_dir"/"$run_name"

mkdir -p "$work_dir"/alignment
mkdir -p "$work_dir"/NanoPlot
mkdir -p "$work_dir"/pycoQC
mkdir -p "$work_dir"/coverage/mosdepth
mkdir -p "$work_dir"/coverage/bedtools
mkdir -p "$work_dir"/SvAnna

#Set up logging
t_stamp=$(date +%H%M_%m%d%y)
exec 3>&1 1>$output_dir/$run_name/$run_name"_"$t_stamp.log 2>&1

echo $(date) >&3
echo "INFO: Output Directory: $work_dir" >&3

## BASECALLING ##

if [ -z "$skip_basecalling" ]
then
echo $(date) >&3
echo "INFO: Basecalling..." >&3

# dorado basecalling with integrated alignmnet (minimap2)
dorado basecaller --device cuda:0 --min-qscore 8 --recursive --reference "$mmi_index" hac@v4.3.0 "$run_dir" > "$work_dir"/alignment/"$run_name".raw.bam

# generate sequencing summary file
echo $(date) >&3
echo "INFO: Generating sequencing summary file..." >&3
dorado summary -v "$work_dir"/alignment/"$run_name".raw.bam > "$work_dir"/alignment/"$run_name".sequencing_summary.tsv

# use samtools to sort, index and generate flagstat file.
# -@ specifies number of threads
else
echo $(date) >&3
echo "INFO: Skipping basecalling" >&3
fi

if [ ! -f "$work_dir"/alignment/"$run_name".bam ]; then
echo "INFO: Sorting and indexing bam file..." >&3
samtools sort \
-@ 20 -o "$work_dir"/alignment/"$run_name".bam "$work_dir"/alignment/"$run_name".raw.bam

#index sorted bam file
samtools index -@ 20 "$work_dir"/alignment/"$run_name".bam

#save stats
echo $(date) >&3
echo "INFO: Generating flagstats..." >&3
samtools flagstat -@ 20"$work_dir"/alignment/"$run_name".bam > "$work_dir"/alignment/"$run_name"_flagstat.txt
fi

# check bam file has been sorted and indexed then clean up unsorted bam file
if [ -s "$work_dir"/alignment/"$run_name".bam ] && [ -s "$work_dir"/alignment/"$run_name".bam.bai]; then
  rm "$work_dir"/alignment/"$run_name".raw.bam
fi

##QC ##

if [ -z "$skip_qc" ]
then
echo $(date) >&3
echo "INFO: Creating summary plots" >&3
  if [ ! -f "$work_dir"/NanoPlot/summary/"$run_name"NanoPlot-report.html ]; then

  echo "INFO: NanoPlot all reads..." >&3
  #plots of run using sequencing summary
  NanoPlot \
  --summary "$work_dir"/alignment/"$run_name".sequencing_summary.tsv \
  --loglength \
  --outdir "$work_dir"/NanoPlot/summary \
  --prefix "$run_name" \
  --threads 20

  #plots of alignment using bam file
  # Use aligned read length not sequence read length
  echo $(date) >&3
  echo "INFO NanoPlot aligned reads..." >&3
  NanoPlot \
  --bam "$work_dir"/alignment/"$run_name".bam \
  --outdir "$work_dir"/NanoPlot/bam \
  --loglength \
  --N50 \
  --prefix "$run_name" \
  --threads 20 \
  --alength 
  fi
fi

##SV CALLING ##
if [ -z $skip_SV ]
then
echo $(date) >&3
echo "INFO: Calling SVs..." >&3
echo "INFO: Phenotype term used - $phenotype_term" >&3

#Sniffles
mkdir "$work_dir"/Sniffles
cd "$work_dir"/Sniffles

sniffles --input "$work_dir"/alignment/"$run_name".bam \
--vcf "$work_dir"/Sniffles/"$run_name"_sniffles.vcf \
--non-germline

# SvAnna annotation prioritising HPO term for myeloid disorders
echo $(date) >&3
echo "INFO: Annotating SVs..." >&3
java -jar ~/Tools/SvAnna/svanna-cli-1.0.4/svanna-cli-1.0.4.jar \
prioritize \
-d /home/matt/Tools/SvAnna/database \
--vcf "$work_dir"/Sniffles/"$run_name"_sniffles.vcf \
--phenotype-term "$phenotype_term" \
--n-threads 10 \
--output-format html \
--prefix "$run_name" \
--report-top-variants 50 \
--out-dir "$work_dir"/SvAnna/

#cd "$pipeline_dir"
fi

##ADAPTIVE SAMPLING ##

# check adaptive sampling output file exists, and get adaptiive sampling data if so
echo $(date) >&3
if [ "$adaptive_sampling" -eq 1 ]
then
  echo "INFO: Adaptive sampling output detected. Processing adaptive sampling data..." >&3
  echo "INFO: Combining adaptive summary files" >&3

#find all adaptive sampling summary files
adaptive_files=$(find "$run_dir" -name 'adaptive_sampling*')

# concatenate adaptive summary files with a single header
awk 'FNR==1 && NR!=1 { while (/^batch_time/) getline; }
    1 {print}' $adaptive_files > "$work_dir"/"$run_name"_combined_adaptive_sampling_summary.csv

#run adaptive sampling analysis script
bash "$pipeline_dir"/SCRIPTS/adaptive.sh -d "$pipeline_dir" \
-n "$run_name" \
-s "$work_dir"/"$run_name"_combined_adaptive_sampling_summary.csv \
-b "$bed_file" \
-w "$work_dir"

# Nanplot on target reads
# use aligned length and filter reads with Q < 8
NanoPlot \
--bam "$work_dir"/alignment/"$run_name"_stop_receiving.sorted.bam \
--outdir "$work_dir"/NanoPlot/on_target \
--loglength \
--N50 \
--prefix "$run_name" \
--threads 20 \
--alength \
--minqual 8 

else 
  echo "INFO: Skipping adaptive sampling analysis." >&3
fi

## COVERAGE ##
#is this necessary?

echo $(date) >&3
echo "INFO: Calculating coverage" >&3

cd "$work_dir"/coverage/mosdepth 

#use mosdepth to calculate depth
mosdepth --by "$bed_file" "$run_name" "$work_dir"/alignment/"$run_name".bam

echo "INFO: Running bedtools" >&3
cd ../bedtools
#get off target reads
#bedtools to find reads in bam file that do and do not not overlap regions in bam
#on target
bedtools intersect -a "$work_dir"/alignment/"$run_name".bam -b "$bed_file" > "$run_name"_on_target.bam
#off target with -v
#might cut this out - doesn't seem that necessary
bedtools intersect -a "$work_dir"/alignment/"$run_name".bam -b "$bed_file" -v > "$run_name"_off_target.bam

echo "INFO: Running samtools" >&3
samtools index "$run_name"_off_target.bam
samtools index "$run_name"_on_target.bam

#get distribution of read lengths
samtools stats "$run_name"_off_target.bam | grep ^RL | cut -f 2- > off_target_len.txt
samtools stats "$run_name"_on_target.bam | grep ^RL | cut -f 2- > on_target_len.txt

#depth and coverage calculations on .tsv output from bedtools
# already done in adaptive.sh??
# Rscript $pipeline_dir/SCRIPTS/coverage_adaptive_panel.r *.tsv $run_name
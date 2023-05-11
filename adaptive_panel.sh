#!/bin/bash

#commands used when processing fastqfile data, following epi2me guides


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
   echo -e "\t-r Path to reference genome mmi index"
   echo -e "\t-a Path to adaptive sampling summary file"
   exit 1 # Exit script after printing help
}

#parse arguments
#TODO change this from getopts
while getopts n:d:b:a:o:r:s:m:q:v:h opt; do
  case "$opt" in
    n) run_name="$OPTARG";;
    d) run_dir="$OPTARG";;
    b) bed_file="$OPTARG";;
    a) adaptive_summary="$OPTARG";;
    o) output_dir="$OPTARG";;
    r) ref_index="$OPTARG";;
    #options to skip steps - useful for isolating/debugging specific aspects
    s) skip_basecalling="$OPTARG";;
    m) skip_alignment="$OPTARG";;
    q) skip_qc="$OPTARG";;
    v) skip_SV="$OPTARG";;
    h) helpFunction;;
  esac
done

#Help if mandatory arguments are empty
if [ -z "$run_name" ] || [ -z "$run_dir" ] || [ -z "$output_dir" ] || \
   [ -z "$ref_index" ] || [ -z "$bed_file" ] || [ -z "$adaptive_summary" ]
then
   echo ""
   echo "ERROR: Missing one or more required arguments. See help message below.";
   helpFunction
fi

#check for running guppyd service before use and exit if running
pid=$( nvidia-smi | grep guppy | awk '{print $5}' )

if [[ "$pid" =~ ^[0-9]+$ ]]; then
  >&2 echo "EXITING: Previous Guppy instance detected. Try: 'sudo service guppyd stop' then retry."
  exit 1
  #kill -9 $pid
 else
   >&2 echo "INFO: No previous Guppy detected, running new analysis..."
fi

#####MAIN######

#create analysis dirs
echo "INFO: Output Directory: $output_dir/$run_name"

mkdir -p "$output_dir"/"$run_name"/alignment
mkdir -p "$output_dir"/"$run_name"/fastq/all
mkdir -p "$output_dir"/"$run_name"/NanoPlot
mkdir -p "$output_dir"/"$run_name"/pycoQC
mkdir -p "$output_dir"/"$run_name"/coverage/mosdepth
mkdir -p "$output_dir"/"$run_name"/coverage/bedtools

#dir variables
pipeline_dir=$(pwd)
work_dir="$output_dir"/"$run_name"


####### BASECALLING ########
#bascall from pod5 files in SUP mode

if [ -z "$skip_basecalling" ] 
then
echo "INFO: Basecalling..."

guppy_basecaller \
--input_path "$run_dir" \
--recursive \
--save_path "$output_dir"/"$run_name"/fastq/all \
--device cuda:0 \
--config dna_r10.4.1_e8.2_400bps_sup.cfg \
--compress_fastq \
--chunks_per_runner 208 \
--gpu_runners_per_device 12 \
--num_callers 4

#merge fastq
cat "$output_dir"/"$run_name"/fastq/all/pass/*.fastq.gz > \
"$output_dir"/"$run_name"/fastq/"$run_name".fastq.gz

fi

####### ALIGNMENT ##########
#alignment step could be combined into guppy command to simplify the pipeline
#maybe more flexible to keep separate?

if [ -z "$skip_alignment" ]
then
echo "INFO: Aligning..."

#align merged fastq to grch38 reference with minimap2
#-a: output SAM file
#-x map-ont: nanopore mode (default)
#using already generated minimap2 indexed reference *.mmi

/opt/ont/guppy/bin/minimap2-2.24 \
-a \
-x map-ont \
"$ref_index" \
"$output_dir"/"$run_name"/fastq/"$run_name".fastq.gz \
 -t 20 > "$output_dir"/"$run_name"/alignment/"$run_name".sam

#use samtools to convert to bam file and sort by position
samtools sort \
-o "$output_dir"/"$run_name"/alignment/"$run_name".bam "$output_dir"/"$run_name"/alignment/"$run_name".sam

#index sorted bam file
samtools index "$output_dir"/"$run_name"/alignment/"$run_name".bam

#save stats
samtools flagstat \
"$output_dir"/"$run_name"/alignment/"$run_name".bam > "$output_dir"/"$run_name"/alignment/"$run_name"_flagstat.txt
fi

######### QC ###############
if [ -z "$skip_qc" ]
then
echo "INFO: Creating summary plots"

pycoQC \
--summary_file "$output_dir"/"$run_name"/fastq/all/sequencing_summary* \
--html_outfile "$output_dir"/"$run_name"/pycoQC/"$run_name"_pycoQC.html \
--bam_file "$output_dir"/"$run_name"/alignment/"$run_name".bam \
--quiet

#plots of run using sequencing summary
NanoPlot \
--summary "$run_dir"/sequencing_summary* \
--loglength \
--N50 \
--outdir "$output_dir"/"$run_name"/NanoPlot/summary \
--prefix "$run_name" \
--threads 20

#plots of alignment using bam file
NanoPlot \
--bam "$output_dir"/"$run_name"/alignment/"$run_name".bam \
--outdir "$output_dir"/"$run_name"/Nanoplot/bam \
--loglength \
--N50 \
--prefix "$run_name" \
--threads 20 \
--alength # Use aligned read length not sequence read length
fi

##### SV CALLING #####
if [ -z $skip_SV ]
then
#cuteSV
mkdir "$work_dir"/CuteSV

#cuteSV for fusion gene detection
#TODO make path to ref genome a variable

cuteSV "$work_dir"/alignment/"$run_name".bam \
~/Tools/ref_genome/grch38/grch38.fa \
"$work_dir"/CuteSV/"$run_name"_cuteSV.vcf \
./ \
--max_cluster_bias_DEL 100 \
--diff_ratio_merging_DEL 0.3 \
--max_size -1 #no uppper limit on SV size 

#Sniffles
mkdir "$work_dir"/Sniffles
cd "$work_dir"/Sniffles

sniffles --input "$work_dir"/alignment/"$run_name".bam \
--vcf "$work_dir"/Sniffles/"$run_name"_sniffles.vcf

cd "$pipeline_dir"
fi

####### ADAPTIVE SAMPLING ##########
#change to specify if adaptive when running command?
#check adaptive sampling output file exists, and get adaptiive sampling data if so

if [ -n "$adaptive_summary" ]
then
  echo "INFO: Adaptive sampling output detected. Processing adaptive sampling data..."

#combine adaptive sampling summary files
#find all adaptive sampling summary files and store in an array
#mapfile -d '' adaptive_files < <(find "$run_dir" -name 'adaptive_sampling*' -print0)
#slice the array to get the header only once in the output file
#{ cat ${adaptive_files[@]:0:1}; grep -v "^batch_time" ${array[@]:1}; } > \
#"$run_name"_combined_adaptive_sampling_summary.csv

#run adaptive sampling analysis script
# TODO try and speed this step up - subsetting bam files takes forever, another way?
  bash "$pipeline_dir"/ALIGNMENT/SCRIPTS/adaptive.sh -d $pipeline_dir \
  -n $run_name \
  -s "$run_name"_combined_adaptive_sampling_summary.csv \
  -b $bed_file \
  -w "$work_dir"
else
  echo "INFO: No adaptive sampling output detected."
fi

###### COVERAGE #####
cd "$work_dir"/coverage/mosepth 

#use mosdepth to calculate depth
mosdepth --by "$bed_file" "$run_name" "$work_dir"/alignment/"$run_name".bam

cd ../bedtools
#get off target reads
#bedtools to find reads in bam file that do and do not not overlap regions in bam
#on target
bedtools intersect -a "$work_dir"/alignment/"$run_name".bam -b "$bed_file" > "$run_name"_on_target.bam
#off target with -v
bedtools intersect -a "$work_dir"/alignment/"$run_name".bam -b "$bed_file" -v > "$run_name"_off_target.bam

samtools index "$run_name"_off_target.bam
samtools index "$run_name"_on_target.bam

#get distribution of read lengths
samtools stats "$run_name"_off_target.bam | grep ^RL | cut -f 2- > off_target_len.txt
samtools stats "$run_name"_on_target.bam | grep ^RL | cut -f 2- > on_target_len.txt



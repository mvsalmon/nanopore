#!/bin/bash

#commands used when processing fastqfile data, following epi2me guides


##TO DO
#ref genome variable?
#logging


#Usage
helpFunction()
{
   echo ""
   echo "Pipeline for analysis of Nanopore runs. Data saved to ~/nanopore_runs/"
   echo ""
   echo "Usage: $0 -h display help -n run_name -d run_dir -o output_dir -b bed_file"
   echo -e "\t-n Name of Nanopore run"
   echo -e "\t-d Path to directory containing Nanopore run data"
   echo -e "\t-b BED file for adaptive sampling analysis" #make this optional
   #TODO echo -e "\t-o Output directory path"
   exit 1 # Exit script after printing help
}

#parse arguments
while getopts n:d:b:a:h opt; do
  case "$opt" in
    n) run_name="$OPTARG";;
    d) run_dir="$OPTARG";;
    b) bed_file="$OPTARG";;
    a) adaptive_summary="$OPTARG";;
    h) helpFunction;;
  esac
done

#Help if mandatory arguments are empty
if [ -z "$run_name" ] || [ -z "$run_dir" ] || [ -z "$bed_file" ]
then
   echo "Some or all of the arguments are empty";
   helpFunction
fi

#check for running guppyd service before use and exit if running
pid=$( nvidia-smi | grep guppy | awk '{print $5}' )

if [[ "$pid" =~ ^[0-9]+$ ]]; then
  >&2 echo "EXITING: Previous Guppy instance detected. Kill and retry"
  exit 1
#   kill -9 $pid
# else
#   >&2 echo "INFO: No previous Guppy detected, running new analysis..."
fi

#####MAIN######

#create analysis dirs
mkdir -p ~/nanopore_runs/"$run_name"/alignment
mkdir -p ~/nanopore_runs/"$run_name"/fastq/all
mkdir -p ~/nanopore_runs/"$run_name"/pycoQC
mkdir -p ~/nanopore_runs/"$run_name"/coverage

#dir variables
pipeline_dir=$(pwd)
work_dir=~/nanopore_runs/"$run_name"

echo "INFO: Working directory $work_dir"
cd $work_dir

#merge fastq files to analysis dirs
cat "$run_dir"/fastq_pass/*.gz > ~/nanopore_runs/"$run_name"/fastq/"$run_name".fastq.gz

####### BASECALLING ########
#bascall from fast5 files in high accuracy mode
#assumes no basecalling during run

echo "INFO: Basecalling..."


#pipe pass and fail fast5 files to guppy
ls "$run_dir"/fast5*/*.fast5 | \
/opt/ont/ont-guppy/bin/guppy_basecaller --save_path "$work_dir"/fastq/all \
--device cuda:0 \
--config dna_r9.4.1_450bps_hac.cfg \
--chunk_size 2000 \
--chunks_per_runner 256 \
--gpu_runners_per_device 2 \
--compress_fastq

# #merge fastq
cat "$work_dir"/fastq/all/pass/*.gz > ~/nanopore_runs/"$run_name"/fastq/"$run_name".fastq.gz

####### ALIGNMENT ##########

#align merged fastq to grch38 reference with minimap2
#-a: output SAM file
#-x map-ont: nanopore mode (default)
#using already generated minimap2 indexed reference *.mmi
echo "INFO: Aligning..."
minimap2 -a -x map-ont \
~/tools/refgenome/seqs_for_alignment_pipelines/grch38/grch38.fasta.gz.mmi \
"$work_dir"/fastq/"$run_name".fastq.gz > "$work_dir"/alignment/"$run_name".sam


#use samtools to convert to bam file and sort by position
samtools sort -o -@ 20 "$work_dir"/alignment/"$run_name".bam "$work_dir"/alignment"$run_name".sam

#index sorted bam file
samtools index "$work_dir"/alignment/"$run_name".bam

#save stats
samtools flagstat "$work_dir"/alignment/"$run_name".bam > "$work_dir"/alignment/"$run_name"_flagstat.txt


######## QC ###############

echo "INFO: Running pycoQC..."

pycoQC \
--summary_file "$work_dir"/fastq/all/sequencing_summary* \
--html_outfile "$work_dir"/pycoQC/"$run_name"_pycoQC.html \
--bam_file "$work_dir"/alignment/"$run_name".bam \
--quiet


####### ADAPTIVE SAMPLING ##########
#change to specify if adaptive when running command?
#check adaptive sampling output file exists, and get adaptiive sampling data if so

if [ -n "$adaptive_summary" ]
then
  echo "INFO: Adaptive sampling output detected. Processing adaptive sampling data..."

#run adaptive sampling analysis script
  bash "$pipeline_dir"/ALIGNMENT/SCRIPTS/adaptive.sh -d $pipeline_dir -n $run_name -s $adaptive_summary -b $bed_file
else
  echo "INFO: No adaptive sampling output detected."
fi


###### COVERAGE #####


mkdir "$work_dir"/coverage/mosdepth
mkdir "$work_dir"/coverage/bedtools

cd "$work_dir"/coverage/mosdepth
#use mosdepth to generate depth
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

cd "$work_dir"
##### CUTESV #####

mkdir ./cuteSV && cd ./cuteSV

#cuteSV for fusion gene detection

cuteSV "$work_dir"/alignment/"$run_name".bam \
~/tools/refgenome/seqs_for_alignment_pipelines/grch38/grch38.fa \
"$run_name".vcf \
./ \
--max_cluster_bias_DEL 100 \
--diff_ratio_merging_DEL 0.3 \
--max_size -1 #no uppper limit on SV size 
#!/bin/bash
#commands used when processing fastqfile data, following epi2me guides
#run from SCRIPTS dir


##TO DO
#pycoQC
#sup basecalling?
#SV  calling

#Usage
helpFunction()
{
   echo ""
   echo "Usage: $0 -h display help -n run_name -d run_dir -o output_dir -b bed_file"
   echo -e "\t-n Name of Nanopre run"
   echo -e "\t-d Directory containng Nanopore run data"
   echo -e "\t-o Directory to write results"
   echo -e "\t-b BED file for adaptive sampling analysis"
   exit 1 # Exit script after printing help
}

#parse arguments
while getopts n:d:o:b:h opt; do
  case "$opt" in
    n) run_name="$OPTARG";;
    d) run_dir="$OPTARG";;
    o) output_dir="$OPTARG";;
    b) bed_file="$OPTARG";;
    h) helpFunction;;
  esac
done

#Help if arguments are empty
if [ -z "$run_name" ] || [ -z "$run_dir" ] || [ -z "$output_dir" ] || [ -z "$bed_file" ]
then
   echo "Some or all of the arguments are empty";
   helpFunction
fi

#####MAIN######

#get pipeline dir
pipeline_dir=$(pwd)

#create analysis dirs
mkdir -p ~/"$output_dir"/"$run_name"/alignment
mkdir -p ~/"$output_dir"/"$run_name"/fastq/pycoQC


########PYCOQC###############
#source conda
source /home/nanopore/miniconda3/etc/profile.d/conda.sh
conda activate pycoQC

pycoQC --summary_file

conda activate py3.8



cd ~/nanopore_runs/"$run_name"

#merge fastqfiles
cat "$run_dir"/fastq_pass/*.gz > ./fastq/"$run_name".fastq.gz

#fastq stats with seqkit
echo "Running seqkit..."
seqkit stats ./fastq/"$run_name".fastq.gz > ./fastq/"$run_name"_seqkit_stats.txt

#align merged fastq to grch38 reference with minimap2
#-a: output SAM file
#-x map-ont: nanopore mode (default)
#using already generated minimap2 indexed reference *.mmi
echo "Aligning..."
minimap2 -a -x map-ont \
~/tools/refgenome/seqs_for_alignment_pipelines/grch38/grch38.fasta.gz.mmi \
./fastq/"$run_name".fastq.gz > ./alignment/"$run_name".sam

cd ./alignment

#use samtools to convert to bam file and sort by position
samtools sort -o "$run_name".bam "$run_name".sam

#index sorted bam file
samtools index "$run_name".bam

#save stats
samtools flagstat "$run_name".bam > "$run_name"_flagstat.txt


cd ../

#check adaptive sampling output file exists, and get adaptiive sampling data if so
adaptive_summary=$(find $run_dir -name adaptive_sampling_*.csv -type f)
if [ -e $adaptive_summary ] ; then

  echo "Adaptive sequencing output detected. Processing adaptive sampling data..."

#run adaptive sampling analysis script
  bash "$pipeline_dir"/SCRIPTS/adaptive.sh -d $pipeline_dir -n $run_name -s $adaptive_summary -b $bed_file

fi

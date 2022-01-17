#!/bin/bash
#commands used when processing fastqfile data, following epi2me guides
#run from SCRIPTS dir


##TO DO
#pycoQC
#sup basecalling?
#SV  calling

run_name=$1
rundir=$2
bed_file=$3

#get pipeline dir
pipeline_dir=$(pwd)

#source conda
source /home/nanopore/miniconda3/etc/profile.d/conda.sh
conda activate py3.8

#create dirs
mkdir -p ~/nanopore_runs/"$run_name"/alignment
mkdir ~/nanopore_runs/"$run_name"/fastq

cd ~/nanopore_runs/"$run_name"

#merge fastqfiles
cat "$rundir"/fastq_pass/*.gz > ./fastq/"$run_name".fastq.gz

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
adaptive_summary=$(find $rundir -name adaptive_sampling_*.csv -type f)
if [ -e $adaptive_summary ] ; then

  echo "Adaptive sequencing output detected. Processing adaptive sampling data..."

#run adaptive sampling analysis script
  bash "$pipeline_dir"/SCRIPTS/adaptive.sh -d $pipeline_dir -n $run_name -s $adaptive_summary -b $bed_file

fi

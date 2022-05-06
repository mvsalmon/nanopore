#!/bin/bash
#commands used when processing fastqfile data, following epi2me guides
#run from SCRIPTS dir


##TO DO
#hac/sup re-basecalling?
#SV  calling

#Usage
helpFunction()
{
   echo ""
   echo "Pipeline for analysis of Nanopore runs. Data saved to ~/nanopore_runs/"
   echo ""
   echo "Usage: $0 -h display help -n run_name -d run_dir -o output_dir -b bed_file"
   echo -e "\t-n Name of Nanopore run"
   echo -e "\t-d Directory containng Nanopore run data"
   echo -e "\t-b BED file for adaptive sampling analysis" #make this optional
   exit 1 # Exit script after printing help
}

#parse arguments
while getopts n:d:b:h opt; do
  case "$opt" in
    n) run_name="$OPTARG";;
    d) run_dir="$OPTARG";;
    b) bed_file="$OPTARG";;
    h) helpFunction;;
  esac
done

#Help if arguments are empty
if [ -z "$run_name" ] || [ -z "$run_dir" ]
then
   echo "Some or all of the arguments are empty";
   helpFunction
fi

#####MAIN######

#create analysis dirs
mkdir -p ~/nanopore_runs/"$run_name"/alignment
mkdir -p ~/nanopore_runs/"$run_name"/fastq/all
mkdir -p ~/nanopore_runs/"$run_name"/pycoQC

#dir variables
pipeline_dir=$(pwd)


#merge fastq files to analysis dirs
#cat "$run_dir"/fastq_pass/*.gz > ~/nanopore_runs/"$run_name"/fastq/"$run_name".fastq.gz

####### BASECALLING ########
#bascall from fast5 files in high accuracy mode
#assumes no basecalling during run

echo "Basecalling..."

/opt/ont/ont-guppy/bin/guppy_basecaller --input_path "$run_dir"/fast5* \
--save_path ~/nanopore_runs/"$run_name"/fastq/all \
--device cuda:0 \
--config dna_r9.4.1_450bps_hac.cfg \
--chunk_size 2000 \
--chunks_per_runner 256 \
--gpu_runners_per_device 2 \
--compress_fastq

#merge fastq
cat ~/nanopore_runs/"$run_name"/fastq/all/fastq_pass/*.gz > ~/nanopore_runs/"$run_name"/fastq/"$run_name".fastq.gz

####### ALIGNMENT ##########
cd ~/nanopore_runs/"$run_name"/

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

######## QC ###############
#source conda
source /home/nanopore/miniconda3/etc/profile.d/conda.sh
conda activate pycoQC

echo "Running pycoQC..."

pycoQC \
--summary_file "$run_dir"/sequencing_summary* \
--html_outfile ~/nanopore_runs/"$run_name"/pycoQC/"$run_name"_pycoQC.html \
--bam_file ~/nanopore_runs/"$run_name"/alignment/"$run_name".bam \
--quiet

conda activate


####### ADAPTIVE SAMPLING ##########
#check adaptive sampling output file exists, and get adaptiive sampling data if so
adaptive_summary=$(find $run_dir -name adaptive_sampling_*.csv -type f)
if [ -e $adaptive_summary ] ; then

  echo "Adaptive sequencing output detected. Processing adaptive sampling data..."

#run adaptive sampling analysis script
  bash "$pipeline_dir"/SCRIPTS/adaptive.sh -d $pipeline_dir -n $run_name -s $adaptive_summary -b $bed_file

fi

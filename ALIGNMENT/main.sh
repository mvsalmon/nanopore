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
if [ -z "$run_name" ] || [ -z "$run_dir" ] || [ -z "$bed_file" ]
then
   echo "Some or all of the arguments are empty";
   helpFunction
fi

#####MAIN######

#get pipeline dir
pipeline_dir=$(pwd)

#create analysis dirs
mkdir -p ~/nanopore_runs/"$run_name"/alignment
mkdir -p ~/nanopore_runs/"$run_name"/fastq/pycoQC

######## QC ###############
#source conda
source /home/nanopore/miniconda3/etc/profile.d/conda.sh
conda activate pycoQC

echo "Running pycoQC..."

pycoQC \
--summary_file "$run_dir"/sequencing_summary* \
--html_outfile ~/"$output_dir"/"$run_name"/fastq/pycoQC/"$run_name"_pycoQC.html \
--quiet

conda activate

####### ALIGNMENT ##########
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

####### ADAPTIVE SAMPLING ##########
#check adaptive sampling output file exists, and get adaptiive sampling data if so
adaptive_summary=$(find $run_dir -name adaptive_sampling_*.csv -type f)
if [ -e $adaptive_summary ] ; then

  echo "Adaptive sequencing output detected. Processing adaptive sampling data..."

#run adaptive sampling analysis script
  bash "$pipeline_dir"/SCRIPTS/adaptive.sh -d $pipeline_dir -n $run_name -s $adaptive_summary -b $bed_file

fi

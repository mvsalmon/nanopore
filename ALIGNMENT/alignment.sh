#! /bin/bash
#commands used when processing fastqfile data, following epi2me guides
#run from SCRIPTS dir

run_name=$1
rundir=$2

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

#extract reads from adaptive sampling region
#/home/nanopore/miniconda3/envs/py3.8/bin/samtools view "$run_name".bam "chr4:51540220-59862419" > "$run_name"_target.sam

cd ../

#check adaptive sampling output file exists, get adaptiive sampling stats if
#make this a separate script file?
adaptive_file=$(find $rundir -name adaptive_sampling_*.csv -type f)
if [ -e $adaptive_file ] ; then

  echo "Adaptive sequencing output detected. Processing adaptive sampling data..."

  mkdir ./adaptive_stats
  cd ./adaptive_stats

  #stats
  conda activate r_env
  Rscript ~/nanopore_runs/SCRIPTS/ALIGNMENT/adaptive_stats.r $adaptive_file $run_name

  conda activate py3.8

  echo "Subseting bam file..."
  python ~/nanopore_runs/SCRIPTS/ALIGNMENT/extract_reads_adaptive.py -b ../alignment/"$run_name".bam -a $adaptive_file -o "$run_name".bam

  #index subsetted files
  for f in *.bam; do
    samtools sort -o $f.sorted.bam $f
    samtools index $f.sorted.bam
  done

  #calculate on target percentages (bedtools)
  #IN PROGRESS
  # for f in *sorted.bam; do
  #   tot_reads=$(samtools flagstat $f | awk '{if(NR==1) print $0}')
  #   AS_reads=$(bedtools coverage -a [BEDFILE] -b $f | awk '{sum+=$5} END{print sum}')

fi

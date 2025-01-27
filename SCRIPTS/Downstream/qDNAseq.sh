#!/bin/bash

# Run QDNAseq on bam files listed in bams.txt

bamfiles="./bams.txt"
work_dir=$(pwd)

for file in bamfiles
do
    name=$basename(${file})
    mkdir ./${name} 
    cd ./${name}
    echo "Processing ${name}..."
    
    Rscript ./QDNAseq.R \
    ${file}

    cd ${work_dir}

done
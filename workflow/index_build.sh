#!/bin/bash -l
#SBATCH -A naiss2026-4-1588
#SBATCH --mem=1760GB
#SBATCH -p memory
#SBATCH --cpus-per-task=128
#SBATCH -t 5-00:00:00
#SBATCH -J NEU_Herps_build
#SBATCH -e /cfs/klemming/projects/snic/sediment_paleogenomics/projects/JAMIE_ALUMBAUGH/outputs/slurm_output/NEU_Herps_build_%A-%a.err
#SBATCH -o /cfs/klemming/projects/snic/sediment_paleogenomics/projects/JAMIE_ALUMBAUGH/outputs/slurm_output/NEU_Herps_build_%A-%a.out

TAXA=NEU_Herps
INPUT_DIR=/cfs/klemming/projects/snic/sediment_paleogenomics/projects/JAMIE_ALUMBAUGH/tools/CPG_BoB/build/new_headers
OUTPUT_DIR=/cfs/klemming/projects/snic/sediment_paleogenomics/projects/JAMIE_ALUMBAUGH/outputs/${TAXA}

mkdir -p $OUTPUT_DIR

#check current versions before running:	
module load bowtie2/2.5.4 
module load samtools/1.20

#if [ ! -f $INPUT_DIR/index_input.csv ]; then
#    touch $INPUT_DIR/index_input.csv
#    ls $INPUT_DIR/*.fna | tr '\n' ',' | sed 's/,$//g'>> $INPUT_DIR/index_input.csv
#   echo "Index input file created."
#    else
#    echo "Index input file already exists."
#fi

#Build the Index
# - references should be FASTA format, file extension can be whatever, e.g. .fa, .fasta, .fna
# - option --large-index is only necessary if your references total over 4 billion bp
bowtie2-build --large-index --threads 64 --seed 260924 $(cat $INPUT_DIR/index_input.csv) $OUTPUT_DIR
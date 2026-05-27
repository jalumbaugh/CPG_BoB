#!/bin/bash -l
#SBATCH -A naiss2025-22-1122
#SBATCH --mem=1760GB
#SBATCH -p memory
#SBATCH --cpus-per-task=128
#SBATCH -t 160:00:00
#SBATCH -J index_verts
#SBATCH -e /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.err
#SBATCH -o /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.out

TAXA=DalarnaPlants
INPUT_DIR=/cfs/klemming/projects/supr/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/database_${TAXA}/masked
OUTPUT_DIR=/cfs/klemming/projects/supr/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/3_bowtie2_${TAXA}/Plants_SE_masked_index

mkdir -p $OUTPUT_DIR

#Run this line to make the bowtie input file if it doesn't already exist:
#ls $INPUT_DIR/*.fna | tr '\n' ',' | sed 's/,$//g'>> $INPUT_DIR/index_input.csv

#check current versions before running:	
module load bowtie2/2.5.4 
module load samtools/1.20

#Build the Index
## Build reference index (called Greenland_GreenlandVerts_index in 4_bowtie2_GreenlandVerts):
# - references should be FASTA format, file extension can be whatever, e.g. .fa, .fasta, .fna
# - option --large-index is only necessary if your references total over 4 billion bp
bowtie2-build --large-index --threads 64 --seed 90210 $(cat $INPUT_DIR/index_input.csv) $OUTPUT_DIR

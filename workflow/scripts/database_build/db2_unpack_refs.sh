#!/bin/bash -l 
#SBATCH -A naiss2025-22-1122
#SBATCH -p shared
#SBATCH -t 02:00:00
#SBATCH -c 32
#SBATCH -J unpack_references
#SBATCH -e /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.err
#SBATCH -o /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.out

TAXA=DalarnaPlants
REF_DIR=database_${TAXA}
OUTPUT_DIR=/cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/$REF_DIR/download

#unpack files from datasets file structure, then rename the downloaded .fna files as TAXONID_GENUSSPECIES.fna

#Use this in place of line 18 if you want to read from a file instead of globbing
#while read file; do
for file in $OUTPUT_DIR/*_refs.zip; do
    # List genome files in zip file; grep out (-o)nly the part matching the (-P)erl regex [^/]*.fna meaning "any number of non-/ characters followed by .fna"
    genomes=($( unzip -l $file ncbi_dataset/data/*/*.fna | grep -oP "[^/]*.fna" ));
    # List .fna files inside the zip (full relative paths)
    fna_files=$(unzip -l "$file" | awk '{print $4}' | grep '\.fna$')
    # Extract only .fna files to OUTPUT_DIR, flattening the directory structure
        for fna in $fna_files; do
        unzip -j "$file" "$fna" -d "$OUTPUT_DIR"
     done
done 
wait
# Add this line if reading from a file instead of globbing
#< $OUTPUT_DIR/downloaded.txt

# Append genus to each extracted file before the .fna
while IFS=$'\t' read -r orig new; do
    for f in "$OUTPUT_DIR"/${orig}_*.fna; do
        # Only rename if file exists
        [ -e "$f" ] && mv -v "$f" "$OUTPUT_DIR/${orig}_${new}.fna"
    done
done < "$OUTPUT_DIR/accession_lookup.tsv"
wait

#print the output FASTA files
ls $OUTPUT_DIR/*.fna > $OUTPUT_DIR/unpacked_refs.txt


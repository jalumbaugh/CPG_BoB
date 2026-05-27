#!/bin/bash -l
#SBATCH -A naiss2025-22-1122
#SBATCH -p shared
#SBATCH -t 20:00:00
#SBATCH -c 32
#SBATCH -J reference_DalarnaPlants
#SBATCH -e /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.err
#SBATCH -o /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.out

TAXA=DalarnaPlants
REF_DIR=database_${TAXA}
OUTPUT_DIR=/cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/$REF_DIR/download

#1. prep taxon_list.txt and/or accession_list.txt
#2. prep accession_lookup.tsv - this is used in the unpacking script, and in steps following this one for masking (GenBank ID, taxon name as genusspecies))

#DO EITHER:
#  datasets download genome taxon --inputfile $OUTPUT_DIR/taxon_list.txt --assembly-level chromosome,complete --include genome --reference --filename $OUTPUT_DIR/taxon_refs.zip
  wait
    #options for --assembly-level include --assembly-level contig,scaffold,chromosome,complete
    #you may remove --reference to download all genome assemblies, but this will download multiple .fnas per taxon! Beware.

#AND/OR
  #download references by GenBank ID (from accession_list.txt)
  #datasets download genome accession --inputfile $OUTPUT_DIR/accession_list.txt --include genome --filename $OUTPUT_DIR/accession_refs.zip
  wait

#LASTLY
  #generate downloaded_references.txt, check to make sure you have all the files you expect
  for file in $(find "$OUTPUT_DIR" -type f -name "*.zip"); 
  do unzip -l $file >> downloaded_references.txt;
  done
  wait 


#!/bin/bash -l

#SBATCH -A naiss2025-22-1122
#SBATCH -p main
#SBATCH --cpus-per-task=64
#SBATCH -t 5:00:00
#SBATCH -J replace_headers
#SBATCH -e /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.err
#SBATCH -o /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.out
#SBATCH --array=0-5

set -euo pipefail

TAXA=DalarnaPlants
INPUT_DIR=/cfs/klemming/projects/supr/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/database_${TAXA}/download
OUTPUT_DIR=/cfs/klemming/projects/supr/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/database_${TAXA}/new_headers
header_files=${OUTPUT_DIR}/newheader_files
input_array=($(<$INPUT_DIR/unpacked_refs.txt))
input_file=${input_array[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTPUT_DIR
mkdir -p $header_files

filename=$(basename "$input_file")            # GCA_000690835.1_Fulmarusglacialis.fna
name=${filename#GCA_}                         # remove leading GCA_ if present
name=${name#GCF_}                             # remove leading GCF_ if present
name=${name%_masked.fna}                      #strip _masked.fna if present
name=${name%.*}                               # strip extension -> 000690835.1_Fulmarusglacialis
IFS=_ read -r acc rest <<< "$name"            # acc=000690835.1  rest=Fulmarusglacialis (rest may contain underscores)
acc_base=${acc%%.*}                           # drop .1 -> 000690835
sample_name="${acc_base}_${rest}_"            # 000690835_Fulmarusglacialis_

#This adds the sample name info to start building a MAP file for db4, need to manually add taxaID from NCBI in second column
echo "${acc_base}_${rest}" >> "$OUTPUT_DIR/acc2taxid_${TAXA}_headers.tsv"

# Extract contig IDs (first token after '>') to contigs_txt
awk '/^>/ { h=substr($0,2); split(h,a,/[^[:alnum:]_.:-]+/); print a[1]; next }' "$input_file" > "$header_files/${sample_name}contigs.txt"
wait

# Replace headers in the fasta with sample_name_contigID (preserves sequence lines)
awk -v s="$sample_name" '
  /^>/ {
    header = substr($0,2)
    # take first token as contig ID (split on whitespace)
    split(header, a, /[[:space:]]+/)
    contig = a[1]
    print ">" s contig
    next
  }
  { print }
' "$input_file" > "$OUTPUT_DIR/${sample_name}newheader.fna"
wait

# Also write a simple list of new headers (sample_contigID lines)
awk -v s="$sample_name" '/^>/ { h=substr($0,2); split(h,a,/ /); print s a[1] }' "$input_file" > "$header_files/${sample_name}new_headers.txt"

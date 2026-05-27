#!/bin/bash -l

#SBATCH -A naiss2025-22-1122
#SBATCH -p main
#SBATCH --cpus-per-task=30
#SBATCH -t 5:00:00
#SBATCH -J compile_acc2taxid
#SBATCH -e /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.err
#SBATCH -o /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.out

set -euo pipefail

TAXA=DalarnaPlants
BASE_DIR=/cfs/klemming/projects/supr/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/database_${TAXA}/new_headers
OUTPUT_DIR=${BASE_DIR}/acc2taxid_output
MAP=$BASE_DIR/acc2taxid_${TAXA}_headers.tsv
DRY_RUN="${DRY_RUN:-0}"
NOMAP="${NOMAP:-NA}"

#OBS! You need to make the MAP file manually from acc2taxid_${TAXA}_headers.tsv from db3 (line 34) for this to work (it pairs the sample names with taxids)

# Usage:
#This script adds taxid info to the acc2taxid headers files created in db3_replaceheaders_refs.sh
#   ./add_taxid_prefix.sh acc2taxid_vert_headers.tsv path/to/*_new_headers.txt
# Optional env:
#   DRY_RUN=1   -> do not overwrite files, just print preview
#   NOMAP=NA    -> value to use when no mapping found

if [[ ! -f "$MAP" ]]; then
  echo "Mapping file not found: $MAP" >&2
  exit 2
fi

# Build map and process each target file found in OUTPUT_DIR
shopt -s nullglob
files=(${BASE_DIR}/newheader_files/*_new_headers.txt)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No matching '*_new_headers.txt' files found in $BASE_DIR/newheader_files" >&2
  shopt -u nullglob
  exit 0
fi
for f in "${files[@]}"; do
  [[ -f "$f" ]] || { echo "Skipping (not a file): $f"; continue; }
  echo "Processing: $f"

  awk -F'\t' -v OFS='\t' -v map="$MAP" -v nomap="$NOMAP" '
    BEGIN {
      while ((getline < map) > 0) {
        if (NF >= 2) {
          mkey = $1
          gsub(/\r/, "", mkey); sub(/[ \t]+$/, "", mkey)
          m[mkey] = $2
        }
      }
      close(map)
    }
    {
  full = $1
  split(full, arr, "_")
  key = arr[1] "_" arr[2]
  gsub(/\r/, "", key); sub(/[ \t]+$/, "", key)
  tax = (key in m ? m[key] : nomap)
  print full, tax
    }
  ' "$f" > "$f.tmp"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "DRY RUN - preview (first 10 lines):"
    head -n 10 "$f.tmp"
    rm -f "$f.tmp"
  else
    cp -p -- "$f" "$f.bak"
    mv -- "$f.tmp" "$f"
    echo "  backed up original to: $f.bak"
  fi
done

# After processing all files, combine them into a single file (all lines from all outputs)
if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "DRY_RUN=1 => skipping combine step"
else
  COMBINED="$OUTPUT_DIR/acc2taxid_${TAXA}.tsv"
  rm -f "$COMBINED"
  shopt -s nullglob
  for g in "$BASE_DIR"/newheader_files/*_new_headers.txt; do
    [[ -f "$g" ]] || continue
    cat "$g" >> "$COMBINED"
  done
  shopt -u nullglob
  if [[ -f "$COMBINED" ]]; then
    echo "Combined files into: $COMBINED ("$(wc -l < "$COMBINED")" lines)"
  else
    echo "No files combined — no matching '*_new_headers.txt' in $BASE_DIR/newheader_files"
  fi
fi
wait

# Duplicate the 1st column (1 & 2 will be identical)
if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "DRY_RUN=1 => skipping prepend of dataset column"
else
  if [[ -f "$COMBINED" ]]; then
    echo "Prepending dataset column to $COMBINED"
    cp -p -- "$COMBINED" "${COMBINED}.orig"
    awk -F"\t" 'BEGIN{OFS="\t"} {print $1, $0}' "$COMBINED" > "${COMBINED}.tmp2" && mv -- "${COMBINED}.tmp2" "$COMBINED"
    echo "Wrote prefixed combined file; backup saved as ${COMBINED}.orig"
  fi
fi
#1.2: downloads fasta files > stage 1 MUST be ran first and separately

import os
import re

#-------------INPUT LISTS-------------
accessions = []
if os.path.exists(config["accession_list"]):
    with open(config["accession_list"], "r") as fhin:
        for line in fhin:
            accessions.append(line.rstrip())

taxa = []
if os.path.exists(config["taxon_list"]):
    with open(config["taxon_list"], "r") as fhin:
        for line in fhin:
            taxa.append(line.rstrip())

valid_taxa = []
for taxon in taxa:
    path = f"build/acc2taxid/accessions/{taxon}.txt"
    if os.path.exists(path):
        with open(path) as fh:
            accession = fh.read().strip()
        if accession:
            valid_taxa.append(taxon)

#-------------WILDCARD CONSTRAINTS-------------
wildcard_constraints:
    taxon="|".join(re.escape(t) for t in valid_taxa) if valid_taxa else "(?!)"

def taxon_zipfile(wildcards):
    return f"build/download/{wildcards.taxon}.zip"

def accession_zipfile(wildcards):
    return f"build/download/{wildcards.accession}.zip"

# -------------LIST OUTPUTS-------------
get_refs_2_outputs = []

get_refs_2_outputs += expand(
    "build/download/{accession}.fna", accession=accessions
)
get_refs_2_outputs += expand(
    "build/download/{taxon}.fna",
    taxon=valid_taxa
)

# -------------RULES-------------
rule download_accessions:
    output:
        fna="build/download/{accession}.fna",
        stamp=temp("logs/user_record/timestamped_{accession}.txt")
    params:
        zipfile=accession_zipfile,
    log:
        "logs/download/{accession}.log"
    shell:
        r"""
        if [ -s {output.fna:q} ]; then
            echo "{output.fna} already exists, skipping download" > {log:q}
            echo "Skipped {wildcards.accession} (reference already downloaded) at $(date)" > {output.stamp:q}
        else
            datasets download genome accession {wildcards.accession} \
                --include genome \
                --filename {params.zipfile:q} > {log:q} 2>&1

            fna_file=$(unzip -Z1 {params.zipfile:q} '*.fna' | head -n1)
            unzip -p {params.zipfile:q} "$fna_file" > {output.fna:q} 2>> {log:q}
            rm -f {params.zipfile:q}

            echo "Downloaded accession {wildcards.accession} at $(date)" > {output.stamp:q}
        fi
        """

rule download_best_acc_for_taxa:
    input:
        accession="build/acc2taxid/accessions/{taxon}.txt"
    output:
        fna="build/download/{taxon}.fna"
    params:
        zipfile=taxon_zipfile,
    log:
        "logs/download_best_acc_for_taxa/{taxon}.log"
    shell:
        r"""
        datasets download genome accession --inputfile {input.accession:q} \
            --include genome \
            --filename {params.zipfile:q} > {log:q} 2>&1

        if [ ! -s {params.zipfile:q} ]; then
            echo "Download did not create {params.zipfile}" >> {log:q}
            exit 1
        fi

        fna_file=$(unzip -Z1 {params.zipfile:q} '*.fna' | head -n1)

        if [ -z "$fna_file" ]; then
            echo "No .fna file found inside {params.zipfile}" >> {log:q}
            rm -f {params.zipfile:q}
            exit 1
        fi

        unzip -p {params.zipfile:q} "$fna_file" > {output.fna:q} 2>> {log:q}
        rm -f {params.zipfile:q}
        """

ruleorder: download_best_acc_for_taxa > download_accessions

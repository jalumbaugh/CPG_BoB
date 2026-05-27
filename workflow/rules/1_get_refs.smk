# 1. Download and unpack references from NCBI from input lists

import os

# -------------LOAD INPUT LISTS-------------
accessions = []
taxa = []

if os.path.exists(config["accession_list"]):
    with open(config["accession_list"], "r") as fhin:
        for line in fhin:
            accessions.append(line.rstrip())

if os.path.exists(config["taxon_list"]):
    with open(config["taxon_list"], "r") as fhin:
        for line in fhin:
            taxa.append(line.rstrip())


# -------------WILDCARD CONSTRAINTS-------------
wildcard_constraints:
    accession="|".join(accessions),
    taxon="|".join(taxa),


# -------------LIST OUTPUTS-------------
get_refs_outputs = []
get_refs_outputs += expand(
    "build/download/{reference}.fna", reference=accessions + taxa
)
get_refs_outputs += expand(
    "logs/user_record/{project_name}_get_references.log",
    project_name=config["project_name"],
)


# --------------RULES-------------
# download from accessions list
rule download_accessions:
    output:
        temp("build/download/{accession}.zip"),
        temp("logs/user_record/timestamped_{accession}.txt"),
    log:
        "logs/download/{accession}.log",
    shell:
        """
        datasets download genome accession {wildcards.accession} \
            --include genome \
            --filename {output[0]} >{log} 2>&1
        echo "Downloaded accession {wildcards.accession} at $(date)" >{output[1]}
        """


# download from taxa list, selecting the best reference genome based on contig N50
rule select_best_acc_from_taxa:
    output:
        temp("build/download/{taxon}.txt"),
    log:
        "logs/select_best_acc_from_taxa/{taxon}.log",
    shell:
        """
        t=$(echo {wildcards.taxon} | tr '_' ' ')
        datasets summary genome taxon "$t" \
            --as-json-lines \
            | dataformat tsv genome --fields accession,assmstats-contig-n50,assminfo-refseq-category \
            | awk -F'\t' 'NR>1 && $2 != "" && tolower($3)=="reference genome" {{print $1 "\t" $2}}' \
            | sort -k2,2nr \
            | head -n1 \
            | cut -f1 >{output} 2>{log}
        echo "Best accession for taxon {wildcards.taxon}: $(cat {output})" >>{log}
        """


# download the best reference genome
rule download_best_acc_for_taxa:
    input:
        "build/download/{taxon}.txt",
    output:
        temp("build/download/{taxon}.zip"),
        temp("logs/user_record/timestamped_{taxon}.txt"),
    log:
        "logs/download_best_acc_for_taxa/{taxon}.log",
    shell:
        """
        datasets download genome accession --inputfile {input} \
            --include genome \
            --filename {output[0]} 2>{log}
        echo "Downloaded best accession for {wildcards.taxon}: $(cat {input}) at $(date)" >{output[1]}
        """


# unpack the fasta files from the downloaded NCBI zip files
rule unpack_references:
    input:
        "build/download/{reference}.zip",
    output:
        "build/download/{reference}.fna",
    log:
        "logs/unpack_references/{reference}.log",
    shell:
        """
        fna_file=$(unzip -Z1 {input} '*.fna' | head -n1)
        unzip -p {input} "$fna_file" >{output} 2>{log}
        """


# make a concatonated log file for this
rule cat_logs:
    input:
        expand("logs/user_record/timestamped_{ref}.txt", ref=accessions + taxa),
    output:
        "logs/user_record/{project_name}_get_references.log",
    shell:
        """
        cat {input} >{output}
        """

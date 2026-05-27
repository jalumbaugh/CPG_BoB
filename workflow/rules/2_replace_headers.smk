# 1. Replace the headers in the raw reference fasta files with the accession name

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
    "build/raw_references/{reference}.fna", reference=accessions + taxa
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
        temp("build/download/{output_accession}.txt"),
    log:
        "logs/user_record/timestamped_{taxon}.log",
    shell:
        """
        t=$(echo {wildcards.taxon} | tr '_' ' ')
        datasets summary genome taxon "$t" \
            --assembly-source refseq \
            --as-json-lines \
            | dataformat tsv genome --fields accession,assmstats-contig-n50,assminfo-refseq-category \
            | awk -F'\t' 'NR>1 && $2 != "" && tolower($3)=="reference genome" {{print $1 "\t" $2}}' \
            | sort -k2,2nr \
            | head -n1 \
            | cut -f1 >{wildcards.output_accession} 2>{log}
        echo "Best accession for taxon {wildcards.taxon}: $(cat {wildcards.output_accession})" >{log}
        """


# download the best reference genome
rule download_best_acc_for_taxa:
    input:
        "build/download/{output_accession}.txt",
    output:
        temp("build/download/{output_accession}.zip"),
        temp("logs/user_record/timestamped_{output_accession}.txt"),
    log:
        "logs/download_best_acc_for_taxa/{output_accession}.log",
    shell:
        """
        datasets download genome accession --inputfile {input} \
            --include genome \
            --filename {output[0]} 2>{log}
        echo "Downloaded accession {wildcards.output_accession}: $(cat {input}) at $(date)" >{output[1]}
        """


# unpack the fasta files from the downloaded NCBI zip files
rule unpack_references:
    input:
        "build/download/{reference}.zip",
    output:
        "build/raw_references/{reference}.fna",
    log:
        "logs/unpack_references/{reference}.log",
    shell:
        """
        fna_file=$(unzip -Z1 {input} '*.fna' | head -n1)
        unzip -p {input} "$fna_file" >{output} 2>{log}
        """


# make a list of downloaded accessions for use in 3_acc2taxid_build.smk
rule accessions_for_acc2taxid:
    output:
        temp("build/acc2taxid/downloaded_accessions.tsv"),
        "build/acc2taxid/acc2taxid_lookup.tsv",
    shell:
        """
        cat {accessions} {output_accession}.txt >{output[0]}
        datasets summary genome accession {output[0]} --as-json-lines \
            | dataformat tsv genome --fields accession,organism-tax-id >{output[1]}
        """


# make a concatonated log file for this project
rule cat_logs:
    input:
        expand(
            "logs/user_record/timestamped_{ref}.txt",
            ref=accessions + taxa + output_accession,
        ),
    output:
        "logs/user_record/{project_name}_get_references.log",
    shell:
        """
        cat {input} >{output}
        """

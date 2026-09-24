#Build your index!

import os
import re

# -------------LOAD INPUT LISTS-------------
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

references = accessions + valid_taxa


# -------------CHOOSE REFERENCE FASTA SOURCE-------------
if config["replace_headers"]:
    ref_fastas = expand("build/new_headers/{reference}.fna", reference=references)
else:
    ref_fastas = expand("build/download/{reference}.fna", reference=references)

# -------------RULES-------------
rule build_bowtie2_index:
    input:
        ref_fastas
    output:
        csv="build/download/index_input.csv",
        expand(
            f"{config['bowtie2_db']}{{ext}}",
            ext=[
                ".1.bt2l",
                ".2.bt2l",
                ".3.bt2l",
                ".4.bt2l",
                ".rev.1.bt2l",
                ".rev.2.bt2l",
            ],
        ),
    params:
        outdir=lambda wildcards: f"build/database_{config['project_name']}"
    log:
        f"{wild}_index_build.log",
    threads: 64
    shell:
        """
        mkdir -p {params.outdir}

        printf "%s\n" {input} | paste -sd, - > {output.csv}

        bowtie2-build --large-index --threads {threads} --seed 90210 $(cat {output.csv}) {params.outdir}
        """


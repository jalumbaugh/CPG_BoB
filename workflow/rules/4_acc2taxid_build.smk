# 3. Build acc2taxid file

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
    path = f"build/.bin/accessions/{taxon}.txt"
    if os.path.exists(path):
        with open(path) as fh:
            accession = fh.read().strip()
        if accession:
            valid_taxa.append(taxon)

references = accessions + valid_taxa


# -------------WILDCARD CONSTRAINTS-------------
wildcard_constraints:
    reference="|".join(re.escape(r) for r in references) if references else "(?!)",


# -------------LIST OUTPUTS-------------
build_acc2taxid_outputs = []
build_acc2taxid_outputs += expand(
    "build/acc2taxid/{project_name}_acc2taxid.tsv.gz",
    project_name=[config["project_name"]],
)
build_acc2taxid_outputs += expand(
    "build/.bin/acc2taxid/{project_name}_acc2taxid_lookup.tsv",
    project_name=[config["project_name"]],
)


# -------------CHOOSE REFERENCE FASTA SOURCE-------------
if config["replace_headers"]:
    ref_fasta_pattern = "build/new_headers/{reference}.fna"
else:
    ref_fasta_pattern = "build/download/{reference}.fna"

# -------------RULES-------------
rule list_ref_headers:
    localrule: 
        True
    input:
        ref_fasta_pattern,
    output:
        temp("build/acc2taxid/{reference}_headers.tsv"),
    log:
        "logs/list_ref_headers/{reference}_headers.log",
    shell:
        """
        awk -v s="{wildcards.reference}" '/^>/ {{h=substr($0,2); split(h,a,/[^[:alnum:]_.:-]+/); print a[1]"\t"s}}' {input} >{output} 2>{log}
        """


# make a list of downloaded accessions for use in 3_acc2taxid_build.smk
rule lookup_for_acc2taxid:
    localrule: 
        True
    input:
        taxon_files=expand("build/.bin/accessions/{taxon}.txt", taxon=taxa),
    output:
        "build/.bin/acc2taxid/{project_name}_acc2taxid_lookup.tsv",
    params:
        all_accessions=lambda wildcards, input: accessions
        + [
            acc.strip()
            for taxon_file in input.taxon_files
            for acc in open(taxon_file)
            if acc.strip()
        ],
    shell:
        """
        datasets summary genome accession {params.all_accessions} --as-json-lines \
            | dataformat tsv genome --fields accession,organism-tax-id,organism-name >{output}
        """


rule add_taxids_accession:
    localrule: 
        True
    input:
        headers="build/acc2taxid/{accession}_headers.tsv",
        lookup=f"build/.bin/acc2taxid/{config['project_name']}_acc2taxid_lookup.tsv",
    output:
        temp("build/acc2taxid/{accession}_acc2taxid.tsv"),
    log:
        "logs/add_taxids/{accession}_acc2taxid.tsv",
    shell:
        """
        awk -F"\t" 'BEGIN{{OFS="\t"}} FNR==NR{{if(FNR>1){{taxid[$1]=$2; name[$1]=$3}}; next}} {{print $1, taxid[$2], name[$2]}}' {input.lookup} {input.headers} >{output} 2>{log}
        """


rule add_taxids_taxon:
    localrule: 
        True
    input:
        headers="build/acc2taxid/{taxon}_headers.tsv",
        taxon_acc="build/.bin/accessions/{taxon}.txt",
        lookup=f"build/.bin/acc2taxid/{config['project_name']}_acc2taxid_lookup.tsv",
    output:
        temp("build/acc2taxid/{taxon}_acc2taxid.tsv"),
    log:
        "logs/add_taxids/{taxon}_acc2taxid.tsv",
    shell:
        """
        acc=$(cat {input.taxon_acc})
        awk -F"\t" -v acc="$acc" 'BEGIN{{OFS="\t"}} FNR==NR{{if(FNR>1){{taxid[$1]=$2; name[$1]=$3}}; next}} {{print $1, taxid[acc], name[acc]}}' {input.lookup} {input.headers} >{output} 2>{log}
        """

rule cat_acc2taxid:
    localrule: 
        True
    input:
        expand("build/acc2taxid/{reference}_acc2taxid.tsv", reference=references),
    output:
        "build/acc2taxid/{project_name}_acc2taxid.tsv",
    log:
        "logs/cat_acc2taxid/{project_name}_acc2taxid_complete.tsv",
    shell:
        """
        cat {input} | sed 's/ /_/g' | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1, $0}}' | awk -F"\t" 'BEGIN{{OFS="\t"}} {{sub(/\..*$/, "", $1); print}}' > {output} 2>{log}
        """

rule zip_acc2taxid:
    localrule: 
        True
    input:
        "build/acc2taxid/{project_name}_acc2taxid.tsv",
    output:
        "build/acc2taxid/{project_name}_acc2taxid.tsv.gz",
    shell:
        """
        gzip -c {input} > {output}
        rm {input}
        """

ruleorder: add_taxids_taxon > add_taxids_accession

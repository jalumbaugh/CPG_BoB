# 3. Build acc2taxid file

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
build_acc2taxid = []
build_acc2taxid += expand(
    "build/acc2taxid/{project_name}_acc2taxid.tsv",
    project_name=config["project_name"],
)


# --------------RULES-------------
# list headers in references, and put in both 1st and 2nd column of acc2taxid file named after reference fna
rule list_ref_headers:
    input:
        "build/download/{reference}.fna",
    output:
        temp("build/acc2taxid/{reference}_headers.tsv"),
    log:
        "logs/list_ref_headers/{reference}_headers.log",
    shell:
        """
        awk -v s="{reference}" '/^>/ { h=substr($0,2); split(h,a,/ /); print s a[1] }' {input} >{output} 2>{log}
        awk -F"\t" 'BEGIN{OFS="\t"} {print $1, $0}' {output}
        """


# Add a third column with the taxids. One rule for accessions, one for taxon named files.
rule add_taxids_accession:
    input:
        "build/acc2taxid/{accession}_headers.tsv",
        "build/acc2taxid/{project_name}_acc2taxid_lookup.txt",
    output:
        temp("build/acc2taxid/{accession}_acc2taxid.tsv"),
    log:
        "logs/add_taxids/{accession}_acc2taxid.tsv",
    shell:
        """
        cp {input[0]} {output}
        awk -F"\t" 'BEGIN{{OFS="\t"}} FNR==NR{{tax[$1]=$2; next}} {{print $1, $2, tax[$1]}}' {input[1]} {input[0]} >{output} 2>{log}
        """


rule add_taxids_taxon:
    input:
        "build/acc2taxid/{taxon}_headers.tsv",
        "build/acc2taxid/{project_name}_acc2taxid_lookup.txt",
    output:
        temp("build/acc2taxid/{taxon}_acc2taxid.tsv"),
    log:
        "logs/add_taxids/{taxon}_acc2taxid.tsv",
    shell:
        """
        cp {input[0]} {output}
        awk -F"\t" 'BEGIN{{OFS="\t"}} FNR==NR{{tax[$3]=$2; next}} {{print $1, $2, tax[$3]}}' {input[1]} {input[0]} >{output} 2>{log}
        """


rule cat_acc2taxid:
    input:
        expand("build/acc2taxid/{reference}_acc2taxid.tsv", reference=accessions + taxa),
    output:
        "build/acc2taxid/{project_name}_acc2taxid.tsv",
    log:
        "logs/cat_acc2taxid/{project_name}_acc2taxid.tsv",
    shell:
        """
        cat {input} >{output} 2>{log}
        """


# this rule needs to:
# copy the contents of {accession}_headers.tsv to {output}, then
# look for the {accession} in the first column of acc2taxid_lookup.txt and add the corresponding 2nd column to the 3rd column of {output}

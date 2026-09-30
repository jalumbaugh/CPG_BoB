# 2. Replace the headers in the raw reference fasta files with the accession name
import os
import re

# -------------INPUT LISTS-------------
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
replace_headers_outputs = []
replace_headers_outputs += expand(
    "build/new_headers/{reference}.fna", reference=references
)
replace_headers_outputs += expand(
    "logs/new_headers/{reference}_new_headers.log", reference=references
)
replace_headers_outputs += expand(
    "build/.bin/new_headers/{reference}_header_lookup.tsv", reference=references
)
replace_headers_outputs += expand(
    "build/reports/{project_name}_header_lookup.tsv",
    project_name=[config["project_name"]],
)

# --------------RULES-------------

rule GENEX_lookup_refs:
    localrule: 
        True
    input:
        fasta="build/download/{reference}.fna",
    output:
        "build/.bin/new_headers/{reference}_header_lookup.tsv",
    shell:
        """
        awk -v s="{wildcards.reference}" '
        /^>/ {{
            header = substr($0,2)
            split(header, a, /[[:space:]]+/)
            contig = a[1]
            print contig "\t" s "_" contig
            next
        }}
        ' "{input.fasta}" >"{output}"
        """

rule cat_GENEX_lookup:
    localrule: 
        True
    input: 
        expand("build/.bin/new_headers/{reference}_header_lookup.tsv", reference=references),
    output:
        "build/reports/{project_name}_header_lookup.tsv",
    shell:
        """
        cat {input} > {output}
        """

rule replace_new_headers:
    localrule: 
        True
    input:
        fasta="build/download/{reference}.fna",
    output:
        fasta="build/new_headers/{reference}.fna",
    log:
        "logs/new_headers/{reference}_new_headers.log",
    shell:
        """
        awk -v s="{wildcards.reference}" '
        /^>/ {{
            header = substr($0,2)
            split(header, a, /[[:space:]]+/)
            contig = a[1]
            print ">" s "_" contig
            next
        }}
        {{ print }}
        ' "{input.fasta}" >"{output.fasta}"
        echo "New headers made for {wildcards.reference} in {output.fasta}" >"{log}"
        """

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
    "build/.bin/download/{reference}_orig.contigs.txt", reference=references
)
replace_headers_outputs += expand(
    "build/new_headers/{reference}.fna", reference=references
)
replace_headers_outputs += expand(
    "logs/orig_headers/{reference}_orig.contigs.log", reference=references
)
replace_headers_outputs += expand(
    "logs/new_headers/{reference}_new_headers.log", reference=references
)


# --------------RULES-------------
rule extract_contigs:
    input:
        "build/download/{reference}.fna",
    output:
        "build/.bin/download/{reference}_orig.contigs.txt",
    log:
        "logs/orig_headers/{reference}_orig.contigs.log",
    shell:
        """
        awk '/^>/ {{ h=substr($0,2); split(h,a,/[^[:alnum:]_.:-]+/); print a[1]; next }}' "{input}" >"{output}"
        echo "Original headers extracted from {wildcards.reference} to {output}" >"{log}"
        """


rule replace_new_headers:
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

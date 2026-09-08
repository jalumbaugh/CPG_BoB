# 2. Replace the headers in the raw reference fasta files with the accession name
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
replace_headers = []
replace_headers += expand("build/download/{reference}.fna", reference=accessions + taxa)
replace_headers += expand(
    "logs/user_record/{project_name}_get_references.log",
    project_name=config["project_name"],
)
replace_headers += expand(
    "build/acc2taxid/{project_name}_acc2taxid_lookup.tsv",
    project_name=config["project_name"],
)


# --------------RULES-------------


# Use acc2taxid_lookup.tsv to find
rule list_ref_headers:
    input:
        "build/download/{reference}.fna",
    output:
        temp("build/acc2taxid/{reference}_headers.tsv"),
    log:
        "logs/list_ref_headers/{reference}_headers.log",
    shell:
        """
        awk -v s="{reference}" '/^>/ { h=substr($0,2); split(h,a,/ /); print s a[1] }' "$input_file" >"$header_files/${sample_name}new_headers.txt"

        datasets download genome accession {wildcards.accession} \
            --include genome \
            --filename {output[0]} >{log} 2>&1
        echo "Downloaded accession {wildcards.accession} at $(date)" >{output[1]}
        """

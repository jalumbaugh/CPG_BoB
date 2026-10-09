# 5. Reports for taxonomy and estimated resources for bowtie2 build

import os
import re

# -------------LIST OUTPUTS-------------
build_preindex_reports_outputs = []
build_preindex_reports_outputs += expand(
    "build/reports/{project_name}_taxonomy_report.tsv",
    project_name=[config["project_name"]],
)
build_preindex_reports_outputs += expand(
    "build/reports/{project_name}_bowtie_build_resources.tsv",
    project_name=[config["project_name"]],
)

# -------------CHOOSE REFERENCE FASTA SOURCE-------------
if config["replace_headers"]:
    ref_fasta_inputs = expand("build/new_headers/{reference}.fna", reference=references)
else:
    ref_fasta_inputs = expand("build/download/{reference}.fna", reference=references)


# -------------RULES-------------
rule taxonomy_report:
    input:
        "build/.bin/acc2taxid/{project_name}_acc2taxid_lookup.tsv",
    output:
        "build/reports/{project_name}_taxonomy_report.tsv",
    localrule: True
    script:
        "gbif_iucn_fetch.py"


rule est_index_resources:
    input:
        ref_fasta_inputs,
    output:
        "build/reports/{project_name}_bowtie_build_resources.tsv",
    localrule: True
    script:
        "est_index_resources.py"

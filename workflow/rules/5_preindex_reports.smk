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
    ref_fasta_dir = "build/new_headers"
else:
    ref_fasta_dir = "build/download"

# -------------RULES-------------
rule taxonomy_report:
    localrule: 
        True
    input:
        "build/.bin/acc2taxid/{project_name}_acc2taxid_lookup.tsv",
    output:
        "build/reports/{project_name}_taxonomy_report.tsv",
    script:
        "gbif_iucn_fetch.py"


rule est_index_resources:
    localrule: 
        True
    input:
        ref_fasta_dir,
    output:
        "build/reports/{project_name}_bowtie_build_resources.tsv",
    script:
        "est_index_resources.py"

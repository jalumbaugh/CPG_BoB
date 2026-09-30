# 1.2: downloads fasta files > stage 1 MUST be ran first and separately

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


# -------------WILDCARD CONSTRAINTS-------------
wildcard_constraints:
    taxon="|".join(re.escape(t) for t in valid_taxa) if valid_taxa else "(?!)",


def taxon_zipfile(wildcards):
    return f"build/download/{wildcards.taxon}.zip"


def accession_zipfile(wildcards):
    return f"build/download/{wildcards.accession}.zip"


# -------------LIST OUTPUTS-------------
download_refs_outputs = []

download_refs_outputs += expand("build/download/{accession}.fna", accession=accessions)
download_refs_outputs += expand("build/download/{taxon}.fna", taxon=valid_taxa)
download_refs_outputs += ["logs/check_duplicate_second_lines/checked.ok"]


# -------------RULES-------------
rule download_accessions:
    localrule: 
        True
    output:
        fna="build/download/{accession}.fna",
        stamp=temp("logs/user_record/timestamped_{accession}.txt"),
    log:
        "logs/download/{accession}.log",
    params:
        zipfile=accession_zipfile,
    shell:
        r"""
        if [ -s {output.fna:q} ]; then
            echo "{output.fna} already exists, skipping download" >{log:q}
            echo "Skipped {wildcards.accession} (reference already downloaded) at $(date)" >{output.stamp:q}
        else
            datasets download genome accession {wildcards.accession} \
                --include genome \
                --filename {params.zipfile:q} >{log:q} 2>&1

            fna_file=$(unzip -Z1 {params.zipfile:q} '*.fna' | head -n1)
            unzip -p {params.zipfile:q} "$fna_file" >{output.fna:q} 2>>{log:q}
            rm -f {params.zipfile:q}

            echo "Downloaded accession {wildcards.accession} at $(date)" >{output.stamp:q}
        fi
        """


rule download_best_acc_for_taxa:
    localrule: 
        True
    input:
        accession="build/.bin/accessions/{taxon}.txt",
    output:
        fna="build/download/{taxon}.fna",
    log:
        "logs/download_best_acc_for_taxa/{taxon}.log",
    params:
        zipfile=taxon_zipfile,
    shell:
        r"""
        datasets download genome accession --inputfile {input.accession:q} \
            --include genome \
            --filename {params.zipfile:q} >{log:q} 2>&1

        if [ ! -s {params.zipfile:q} ]; then
            echo "Download did not create {params.zipfile}" >>{log:q}
            exit 1
        fi

        fna_file=$(unzip -Z1 {params.zipfile:q} '*.fna' | head -n1)

        if [ -z "$fna_file" ]; then
            echo "No .fna file found inside {params.zipfile}" >>{log:q}
            rm -f {params.zipfile:q}
            exit 1
        fi

        unzip -p {params.zipfile:q} "$fna_file" >{output.fna:q} 2>>{log:q}
        rm -f {params.zipfile:q}
        """


ruleorder: download_best_acc_for_taxa > download_accessions


rule check_duplicate_second_lines:
    localrule: 
        True
    input:
        expand("build/download/{accession}.fna", accession=accessions),
        expand("build/download/{taxon}.fna", taxon=valid_taxa),
    output:
        touch("logs/check_duplicate_second_lines/checked.ok"),
    log:
        "logs/check_duplicate_second_lines/warnings.txt",
    run:
        import os
        import glob

        os.makedirs("logs/check_duplicate_second_lines", exist_ok=True)
        seen = {}
        duplicates = []
        with open(log[0], "w") as fout:
            for f in input:
                if not os.path.exists(f) or os.path.getsize(f) == 0:
                    continue
                with open(f) as fin:
                    fin.readline()
                    second_line = fin.readline().rstrip("\n")
                if second_line in seen:
                    duplicates.append((f, seen[second_line], second_line))
                    fout.write(
                        f"ERROR: {f} and {seen[second_line]} share the same second line: {second_line}\n"
                    )
                else:
                    seen[second_line] = f
            if duplicates:
                fout.write("Removing build files...\n")
                for path in glob.glob("build/.bin/accessions/*"):
                    if os.path.isfile(path):
                        os.remove(path)
                for path in glob.glob("build/.bin/acc2taxid/*"):
                    if os.path.isfile(path):
                        os.remove(path)
                for path in glob.glob("build/acc2taxid/*_acc2taxid.tsv"):
                    if os.path.isfile(path):
                        os.remove(path)
                for path in glob.glob("build/acc2taxid/*_headers.tsv"):
                    if os.path.isfile(path):
                        os.remove(path)
        if duplicates:
            raise ValueError(
                "DUPLICATE FILES DETECTED! See logs/check_duplicate_second_lines/warnings.txt"
            )
        with open(output[0], "w") as ok:
            ok.write("ok\n")

# 1.1: Download and unpack references from NCBI from input lists

import os
import re

# -------------INPUT LISTS-------------
taxa = []
if os.path.exists(config["taxon_list"]):
    with open(config["taxon_list"], "r") as fhin:
        for line in fhin:
            taxa.append(line.rstrip())


# -------------WILDCARD CONSTRAINTS-------------
wildcard_constraints:
    taxon="|".join(re.escape(t) for t in taxa) if taxa else "(?!)",


# -------------LIST OUTPUTS-------------
acc_from_taxa_outputs = []
acc_from_taxa_outputs += expand("build/.bin/accessions/{taxon}.txt", taxon=taxa)
acc_from_taxa_outputs += expand(
    "build/acc2taxid/{project_name}_accfromtax_report.tsv",
    project_name=[config["project_name"]],
)


# -------------RULES-------------
rule select_best_acc_from_taxa:
    output:
        "build/.bin/accessions/{taxon}.txt",
    log:
        "logs/select_best_acc_from_taxa/{taxon}.log",
    shell:
        r"""
        set +o pipefail
        t=$(echo {wildcards.taxon:q} | tr '_' ' ')

        datasets summary genome taxon "$t" \
            --as-json-lines \
            | dataformat tsv genome --fields accession,assmstats-contig-n50,assminfo-refseq-category \
                >{log:q}

        awk -F'\t' 'NR>1 && $2 != "" && tolower($3)=="reference genome" {{print $1 "\t" $2}}' {log:q} \
            | sort -k2,2nr \
            | head -n1 \
            | cut -f1 >{output:q}

        echo "Taxon queried: $t" >>{log:q}
        echo "Selected accession: $(cat {output:q})" >>{log:q}
        """


rule acc_from_taxa_report:
    input:
        expand("build/.bin/accessions/{taxon}.txt", taxon=taxa),
    output:
        "build/acc2taxid/{project_name}_accfromtax_report.tsv",
    log:
        "logs/acc_from_tax/{project_name}_warnings.txt",
    shell:
        r"""
        >{output:q}
        >{log:q}

        for i in {input:q}; do
            taxon=$(basename "$i" .txt)

            if [[ -s "$i" ]]; then
                echo "$taxon= $(cat "$i")" >>{output:q}
            else
                echo "$taxon= no NCBI accession found" >>{output:q}
            fi
        done

        sort -k2,2 {output:q} -o {output:q}

        for i in {input:q}; do
            for j in {input:q}; do
                if [[ "$i" < "$j" ]]; then
                    if [[ -s "$i" && -s "$j" ]]; then
                        if cmp -s "$i" "$j"; then
                            warning="WARNING! These files share the same accession: $i $j"
                            echo "$warning" >&2
                            echo "$warning" >>{log:q}
                            echo "Removing build files..."
                            rm -f build/.bin/accessions/*
                            exit 1
                        fi
                    fi
                fi
            done
        done
        """

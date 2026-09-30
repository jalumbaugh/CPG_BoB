# Build your index!
import os
import re

# -------------CHOOSE REFERENCE FASTA SOURCE-------------
if config["replace_headers"]:
    ref_dir = "build/new_headers"
else:
    ref_dir = "build/download"

# -------------LIST OUTPUTS-------------
build_bowtie2_index_outputs = []
build_bowtie2_index_outputs += ["build/index/index_input.csv"]
build_bowtie2_index_outputs += [
    directory(config["build_bowtie2_index"]["output_dir"] + "/" + config["project_name"])
]
build_bowtie2_index_outputs += [
    f"logs/build_bowtie2_index/{config['project_name']}.log"
]


# -------------RULES-------------
rule make_input_list:
    localrule: 
        True
    input:
        ref_dir,
    output:
        "{ref_dir}/index_input.csv",
    shell:
        r"""
        if [ ! -f {output:q} ]; then
            touch {output:q}
            readlink -f {input}/*.fna | tr '\n' ',' | sed 's/,$//g' >>{output:q}
            echo "Index input file created."
        else
            echo "Index input file already exists."
        fi
        """


rule build_bowtie2_index:
    input:
        "build/index/index_input.csv",
    output:
        outdir=directory(
            config["build_bowtie2_index"]["output_dir"] + "/" + config["project_name"]
        ),
    log:
        f"logs/build_bowtie2_index/{config['project_name']}.log",
    params:
        seed=config["build_bowtie2_index"]["seed"],
        bowtie2_threads=config["build_bowtie2_index"]["bowtie2_threads"],
    shell:
        r"""
        mkdir -p "$(dirname {log})"
        mkdir -p {output.outdir}

        {{
            module load bowtie2/2.5.4
            module load samtools/1.20

            bowtie2-build \
                --large-index \
                --threads {params.bowtie2_threads} \
                --seed {params.seed} \
                "$(cat {input})" \
                {output.outdir}
        }} > {log} 2>&1 || {{
            rm -f build/index/*.sa build/index/*.tmp
            exit 1
        }}
        """







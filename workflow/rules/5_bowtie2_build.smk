# Build your index!

import os
import re

# -------------LIST OUTPUTS-------------
build_bowtie2_index_outputs = []
build_bowtie2_index_outputs += expand("{ref_dir}/index_input.csv")
build_bowtie2_index_outputs += expand(
    outdir=directory(
        config["build_bowtie2_index"]["output_dir"] + "/" + config["project_name"]
    ),
)

# -------------CHOOSE REFERENCE FASTA SOURCE-------------
if config["replace_headers"]:
    ref_dir = expand("build/new_headers")
else:
    ref_dir = expand("build/download")


# -------------RULES-------------
rule make_input_list:
    input:
        ref_dir,
    output:
        "{ref_dir}/index_input.csv",
    shell:
        r"""
        if [ ! -f {output:q} ]; then
            touch {output:q}
            ls {input}/*.fna | tr '\n' ',' | sed 's/,$//g' >>{output:q}
            echo "Index input file created."
        else
            echo "Index input file already exists."
        fi
        """


rule build_bowtie2_index:
    input:
        "{ref_dir}/index_input.csv",
    output:
        outdir=directory(
            config["build_bowtie2_index"]["output_dir"] + "/" + config["project_name"]
        ),
    params:
        seed=config["build_bowtie2_index"]["seed"],
        bowtie2_threads=config["build_bowtie2_index"]["bowtie2_threads"],
    shell:
        r"""
        mkdir -p {output.outdir}

        module load bowtie2/2.5.4
        module load samtools/1.20

        bowtie2-build \
            --large-index \
            --threads {params.bowtie2_threads} \
            --seed {params.seed} \
            $(cat {input}) \
            {output.outdir}
        """

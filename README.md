# **CPG_BoB**
Centre for Palaeogenetics Builder for Bowtie2: 
A bowtie2 index building companion in SnakeMake for large metagenomic databases

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

## Setup: Before starting the workflow
<b>Set up the environment for BoB:</b>
* from the CPG_BoB main directory, run `conda env create -f workflow/envs/BoB_env.yaml`
* then activate the environment `conda activate BoB_env`

<b>Build your taxon and/or accession input files</b>
* you can add any taxon you are interested to your taxon list, as BoB will not fail if no accession exists.
* the taxon list should be formatatted like this, with underscores between the genus and species names:
```
Lepidocyrtus_curvicollis
Drosophila_melanogaster
Sophophora_melanogaster
Hyla_arborea
```
* the accession list should be formatatted like this.
```
GCA_036971115.1
GCA_036687305.1
```
* examples of the input lists can be found under `.test`.

<b>Complete the config file</b>
* The default config is in `config/config.yml`.
* Put the list paths into your config file under PIPELINE RESOURCES.
* Don't forget to change the `project_name` to whatever you want your index to be called.
* Select the rules to use.

## Running the pipeline
* To run the pipeline on a local machine, go to the base directory for CPG_BoB and run `snakemake -s workflow/Snakefile --configfile config/config.yaml -c 1 -p`
    * `--configfile` explicitly states the config file to use. If you are using the default config in `config/config.yml`, you can remove this flag.
    *  `-c` sets the number of cores to use (1 is usually enough)
    *  `-p` prints the output to the screen
    *  `-n` can be added to the end of the run command for a dry run that does not produce output. This is highly recommended!

* To run on a cluster, do the following:
    *  `salloc -A [projectID] -p shared -t 1-00:00:00 --cpus-per-task=8`
    *  `ssh [assigned_node]`
    *  `ml tmux`
    *  `tmux new-session -s [session_name]`
    *  `cd /pathto/CPG_BoB`
    *  `conda activate BoB_env`
    *  `snakemake -s workflow/Snakefile --configfile config/config.yaml -c 1 -p`

* If you are using a taxon list (or a taxon AND accession list) as input:
  * You must first run `get_references_stage1= True` with all the other rules as `= False` before doing anything else.
  * Once you have ran  `get_references_stage1`, check the output `[project_name]_accfromtax_report.tsv` to see which taxa succeeded or failed, e.g.:
      ```
      Lepidocyrtus_curvicollis= GCA_964276635.1
      Drosophila_melanogaster= GCA_000001215.4
      Sophophora_melanogaster= GCA_000001215.4
      Hyla_arborea= no NCBI accession found
      WARNING! These files share the same accession: build/acc2taxid/accessions/Drosophila_melanogaster.txt build/acc2taxid/accessions/Sophophora_melanogaster.txt
      ```
    * If you see something like the WARNING above, then you have accidentally added multiple taxon names that are synonyms or otherwise have the same NCBI accession. You will need to remove one of these from your input taxon list before continuing.
    * If you have generated a blank `.txt`, then you have simply added a blank line to the end of your input taxon list.
      
* If you are only using and accession list as input, you may proceed directly to the full pipeline (i.e. `get_references_stage2` and onwards).

## Development Plans
   * Completion of 4_bowtie2_build with parameters for cluster submission
   * parameter options for seleting the best accession for taxa (currently defaults to reference genomes or those with the highest contig N50):
      * highest assmstats-contig-l50
      * most recent assminfo-release-date
      * refseq only

### Author
Jamie Alumbaugh <br>
Centre for Palaeogenetics, Stockholm
 
### Ackowledgements
Special thanks to NBIS (National Bioinformatics Infrastructure Sweden) for their Snakemake BYOC (bring-your-own-code) Workshop of Spring 2026. 

The basic github framework of this pipeline is based on [this template](https://github.com/snakemake-workflows/snakemake-workflow-template).

This pipeline was also developed with documentation search and coding assistance from [Snakemake Guru AI](https://snakemake.readthedocs.io/en/stable/#).
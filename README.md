# **CPG_BoB**
Centre for Palaeogenetics Builder for [Bowtie2](https://github.com/BenLangmead/bowtie2/blob/master/README.md)]: 
A bowtie2 index building companion in SnakeMake for large metagenomic databases

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

## Setup: Before starting the workflow
<b>Set up the environment for BoB:</b>
* from the CPG_BoB main directory, run `conda env create -f workflow/envs/BoB_env.yaml`
* then activate the environment `conda activate BoB_env`

<b>Build your taxon and/or accession input files</b>
* you can add any taxon you are interested to your taxon list, as BoB will not fail if no accession exists.
* the list may include species, genera, families, or higher taxonomic orders.
* the taxon list should be formatatted like this, with underscores between the genus and species names:
```
Lepidocyrtus_curvicollis
Drosophila_melanogaster
Sophophora_melanogaster
Hyla_arborea
Hyla
```
* the accession list should be formatatted like this.
```
GCA_036971115.1
GCA_036687305.1
```
* examples of the input lists can be found under `.test`.

<b>Complete the config file</b>
* The default config is in `config/config.yml`.
* Change the `project_name` to whatever you want your index to be called.
* Put the list paths into your config file under PIPELINE RESOURCES.
* Mark the rules to want to use as `True`:
    * `get_accessions_from_taxa` searches NCBI for the GenBank reference genome with the highest contig N50 value for each name in your `taxon_list` and pulls its accession number. It will generate a report of the results.
    * `download_references` downloads and unpacks the files matching the accession IDs provided in your `accession_list` and in the .bin/accessions folder produced by `get_accessions_from_taxa`. This rule also creates a lookup file that lists the taxa included in the download folder, as well as their NCBI accession and taxid values.
    * `replace_headers` will replace the headers in the downloaded fasta files with the pattern "taxon/reference_contigID".
        *This rule is optional. If you do not use it, BoB will use the output from `download_references` to build your index.
        *If you do use `replace_headers`, it must remain marked as `True` in the config when running the following steps.
        *Replacing headers is helpful to prevent terminal errors in places where contigIDs may be reused, e.g. when mapping against multiple indexes that may potentially include the same reference fasta multiple times. 
    * `build_acc2taxid` will create your acc2taxid file, which is needed for programs like ngsLCA to correlate mapped contigs with taxonomy. It will also produce a report that includes the accession ID, taxid, full taxonomic path, and global distribution of each organism in your database based on [GBIF](https://www.gbif.org/).
    * `build_bowtie_index` is in development!

## Running the pipeline
* To run the pipeline on a local machine, go to the base directory for CPG_BoB and run `snakemake -s workflow/Snakefile --configfile config/config.yaml --directory ../test -c 1 -p`
    * `--configfile` explicitly states the config file to use. If you are using the default config in `config/config.yml`, you can remove this flag.
    * `--directory` is an optional flag that allows you to put the database build elsewhere instead of in the default location (i.e. CPG_BoB/build).
    *  `-c` sets the number of cores to use (1 is usually enough)
    *  `-p` prints the output to the screen
    *  `-n` can be added to the end of the run command for a dry run that does not produce output. This is highly recommended!

* To run on a cluster, do the following:
    *  `salloc -A [projectID] -p shared -t 1-00:00:00 --cpus-per-task=3`
    *  `ssh [assigned_node]`
    *  `ml tmux`
    *  `tmux new-session -s [session_name]`
    *  `cd /pathto/CPG_BoB`
    *  `conda activate BoB_env`
    *  `snakemake -s workflow/Snakefile --configfile config/config.yaml -c 1 -p`

* If you are using a taxon list (or a taxon AND accession list) as input:
  * Again: you must first run `get_accessions_from_taxa= True` with all the other rules as `= False` before doing anything else.
  * Once you have ran `get_accessions_from_taxa`, check the output `[project_name]_accfromtax_report.tsv` to see which taxa succeeded or failed, e.g.:
      ```
      Lepidocyrtus_curvicollis= GCA_964276635.1
      Drosophila_melanogaster= GCA_000001215.4
      Sophophora_melanogaster= GCA_000001215.4
      Hyla_arborea= no NCBI accession found
      ```
    * For taxa where no NCBI accession was found, you may want to look for a genus or family level representative if there are no others currently in your list. Simply update the taxon list with the genus or family names desired and rerun `get_references_stage1= True`. There is no need to remove the taxa that did not return accessions.
    > [!WARNING]
    > If this step fails, you may have accidentally included two taxa which are synonymous or otherwise have the same NCBI accession. You will need to remove all but one of these from your input lists before continuing. The warning report can be found in `logs/acc_from_tax/{project_name}_warnings.txt`.

    > If you have generated a blank `.txt`, then you have simply added a blank line to the end of your input taxon list.
      
    * If you are only using and accession list as input, you may proceed directly to the full pipeline (i.e. `download_references` and onwards).
        > [!NOTE]
        > If the pipeline hits an unexplained error partway through downloading, try to just run it again: this has occured before when .zip files are very large and internet connections are unstable. 
        > [!WARNING]
        > There is an additional check built into `download_references` that will kill the pipeline and produce the error file `logs/check_duplicate_second_lines/warnings.txt` if mulitple downloaded reference fastas include identical sequences in the first contig. This is to help prevent users from indexing the same reference twice (e.g. a CGA and GCF from the same accession ID). To fix your build: remove all but one of the duplicated .fna files (you may leave the others in place), remove the duplicate entries from your input lists, and start over from the top of the pipeline.


## Development Plans
   * Completion of 5_bowtie2_build with parameters for cluster submission

### Author
Jamie Alumbaugh <br>
Centre for Palaeogenetics, Stockholm
 
### Ackowledgements
Special thanks to NBIS (National Bioinformatics Infrastructure Sweden) for their Snakemake BYOC (bring-your-own-code) Workshop of Spring 2026. 

The basic github framework of this pipeline is based on [this template](https://github.com/snakemake-workflows/snakemake-workflow-template).

This pipeline was also developed with documentation search and coding assistance from [Snakemake Guru AI](https://snakemake.readthedocs.io/en/stable/#).
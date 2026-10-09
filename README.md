# **CPG_BoB**
Centre for Palaeogenetics Builder for [Bowtie2](https://github.com/BenLangmead/bowtie2/blob/master/README.md): 
an index building companion in SnakeMake for metagenomic databases

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

## Setup: Before starting the workflow
<b>Set up the environment for BoB:</b>
* from the CPG_BoB main directory, run `conda env create -f workflow/envs/BoB_env.yaml`
* then activate the environment `conda activate BoB_env`

<b>Build your taxon and/or accession input files</b>
* you can add any taxon you are interested to your taxon list, as BoB will not fail if no accession exists.
* the list may include species, genera, families, or higher taxonomic orders.
* the taxon list should be formatted like this, with underscores between the genus and species names:
```
Drosophila_melanogaster
Clunio_africanus
Clunio
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

    * `get_accessions_from_taxa` searches NCBI for the GenBank reference genome with the highest contig N50 value for each name in your `taxon_list` and pulls its accession number. It will generate a report of the results called `[project_name]_accfromtax_report.tsv` that tells you which accession ID was selected for each taxon.

    * `download_references` downloads and unpacks the files matching the accession IDs provided in your `accession_list` and in the .bin/accessions folder produced by `get_accessions_from_taxa`. 

    * `replace_headers` will replace the headers in the downloaded fasta files with the pattern "{reference}_contigID".
        * <b>This rule is optional. If you do not use it, BoB will use the output from `download_references` to build your index.</b>
        * If you do use `replace_headers`, it must remain marked as `True` in the config when running all steps that follow it.
        * Replacing headers is helpful to prevent terminal mapping errors in places where contigIDs may be reused, e.g. when mapping against multiple indexes that may potentially include the same reference fasta multiple times. It is also useful when the header lines themselves are too long, causing the [header section of the post-mapping .sam files to exceed 2GB and become unable to be coerced into a .bam file](https://github.com/samtools/htslib/issues/1420). The same problem occurs when the headers are short but are many... in which case you will need another solution (forthcoming in this pipeline).
        * This part of the pipeline will also generate a report called `[project_name]_header_lookup.tsv`, where the first column is the original header name, and the second is the new header. This file is very useful to keep around, particularly if you have generated .bed files (such as with the [GENEX pipeline](https://github.com/NikolayOskolkov/MCWorkflow) or [metaJAM](https://github.com/NathanACO/metaJAM)) based on the original header names and need to swap them out for the new headers.

    * `build_acc2taxid` will create and compress your acc2taxid file, which is needed for programs like ngsLCA to correlate mapped contigs with taxonomy. 
    
    * `build-preindex_reports` will generate two reports:
        * a report that includes the accession ID, taxid, full taxonomic path, and global distribution of each organism in your database based on [GBIF](https://www.gbif.org/) called `[project_name]_taxonomy_report.tsv`
        * in beta: a report called `[project_name]_bowtie_build_resources.tsv` that will estimate the computational resources needed for bowtie2 build, which you can set in the snakemake profile (profiles/config.yaml). It will also estimate the size of the final index--- make sure you have enough storage space for all the files!

    * `build_bowtie_index` will build your index in the output folder specified. 

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

* The only rule which does not run locally on the interactive node is `build_bowtie2_index`. Before you run this:
    * use `[project_name]_bowtie_build_resources.tsv` to help you adjust the SLURM settings in `profiles/default/config.yaml`. 
    * This feature is still being tested, so adjust your settings with that in mind.

### If you are using a taxon list (or a taxon AND accession list) as input:
  * While the pipeline will run with everything set to `True`, it is best to first run `get_accessions_from_taxa= True` with all the other rules as `= False`.
    * This way you can modify the input lists to account for any taxa that were not found on NCBI before proceeding to the download stage. 
  * Once you have ran `get_accessions_from_taxa`, check the output `[project_name]_accfromtax_report.tsv` to see which taxa succeeded or failed, e.g.:
      ```
      Drosophila_melanogaster= GCA_000001215.4
      Clunio= GCA_900005825.1
      Clunio_africanus= no NCBI accession found
      ```
    * For taxa where no NCBI accession was found, you may want to look for a genus or family level representative if there are no others currently in your list. Simply update the taxon list with the genus or family names desired and rerun `get_references_stage1= True`. There is no need to remove the taxa that did not return accessions.

> [!NOTE]
> If you have generated a blank `.txt` files, then you have simply added a blank line to the end of your input taxon list.

> [!WARNING]
> If this step fails, you may have accidentally included two taxa which are synonymous or otherwise have the same NCBI accession. You will need to remove all but one of these from your input lists before continuing. The warning report can be found in `logs/acc_from_tax/{project_name}_warnings.txt`.


### If you are only using an accession list as input, or have completed `get_accessions_from_taxa`
  * Set `download_references` and all other rules you wish to run to `True`.
   
> [!NOTE]
> If the pipeline hits an unexplained error partway through downloading, try to run it again: this has occured before when .zip files are very large and internet connections are unstable.

> [!WARNING]
> There is an additional check built into `download_references` that will kill the pipeline and produce the error file `logs/check_duplicate_second_lines/warnings.txt` if mulitple downloaded reference fastas include identical sequences in the first contig. This is to help prevent users from indexing the same reference twice (e.g. a CGA and GCF from the same accession ID). <br> To fix your build: remove all but one of the duplicated .fna files (you may leave the others in place), remove the duplicate entries from your input lists, and start over from the top of the pipeline.

  * After the download stage is complete, `replace_headers` will replace the contig names in each fasta file with the pattern described above if it was set to `True`. The report `[project_name]_header_lookup.tsv` will look something like this, with the original header followed by the new version:
      ```
      CP146036.1	GCA_036971115.1_CP146036.1
      CP146037.1	GCA_036971115.1_CP146037.1
      CP146038.1	GCA_036971115.1_CP146038.1
      (continued)
      ```

  * For `build_acc2taxid`, the acc2taxid file (`[project_name]_acc2taxid.tsv.gz`) will be compressed and formatted in 4 columns to be compatible with other indexes prepared for use in [metaJAM](https://github.com/NathanACO/metaJAM). Here are some example lines from the file made from the .test input lists: 
      ```
      GCA_036971115	GCA_036971115.1_CP146036.1	90675	Camelina_sativa
      GCA_036687305	GCA_036687305.1_CP136102.1	3656	Cucumis_melo
      Drosophila_melanogaster_AE014298	Drosophila_melanogaster_AE014298.5	7227	Drosophila_melanogaster
      Clunio_CVRI01006194	Clunio_CVRI01006194.1	568069	Clunio_marinus
      (continued)
      ```

  * The rule `build_preindex_reports` will generate the report `[project_name]_taxonomy_report.tsv`:
      ```Assembly Accession	Organism Taxonomic ID	class	order	family	genus	species	countries	country_source
      GCA_036971115.1	90675	Magnoliopsida	Brassicales	Brassicaceae	Camelina	sativa	Austria; Belgium; Switzerland (cont.)	gbif
      GCA_900005825.1	568069	Insecta	Diptera	Chironomidae	Clunio	marinus	Belgium; Germany; Denmark  (cont.)	gbif
      GCA_036687305.1	3656	Magnoliopsida	Cucurbitales	Cucurbitaceae	Cucumis	melo	United Arab Emirates; Argentina; Austria (cont.)	gbif
      GCA_000001215.4	7227	Insecta	Diptera	Drosophilidae	Drosophila	melanogaster	Argentina; American Samoa; Austria (cont.)	gbif
      (continued)
      ```

  * ... as well as `[project_name]_bowtie_build_resources.tsv`:
      ```
      Folder=	build/new_headers
      Number of files=	4
      Input .fna size (GiB)=	1.20
      Uncompressed .fna size (GiB)=	1.20
      Number of contigs=	25665

      Time Est. model=	0.0552 * GiB^1.24 h
      Core Est. model=	-0.044*GiB + 32.6 cores
      Index Est. model=	1.89 GiB/GiB input + 0.119 GiB/M contigs

      SNAKEMAKE SETTINGS (profiles/config.yaml)
      Suggested mem_mb=	15116.0
      Suggested runtime (minutes)=	12.46
      Suggested cpus_per_task=	33

      SLURM SETTINGS EQUIVALENT
      Suggested --mem (Mb)=	15117
      Suggested --time=	0-00:12:28
      Suggested --cpus-per-task=	33
      NOTE: there is no need to set the number of cpus if your HPC grants an entire node on a memory partition.

      Est. Index Size (GiB)=	2.3
      NOTE: Check free disk space at the output path before running bowtie2 build.

	  Note: 1 GiB is outside the calibrated range (68-152 GiB); the estimate is an extrapolation.; [cpus/NOTE] there is no need to set the number of cpus if your HPC grants an entire node on a memory partition.; [index/INFO] Check free disk space at the output path before running bowtie2-build --large-index.
      ```
      
  * Once you are happy with your outputs, enable `build_bowtie2_index` and set the output path. Remember to modify the settings under `profiles/config.yaml` before running the pipeline again.


## Development Plans
    * More parameters: 
        * user settings to adjust how accession IDs are selected for taxon names
        * more options for how to re-name replaced headers
    * Stand-alone (simple) tools:
        * .bed file header conversion tool that uses `[project_name]_header_lookup.tsv` to swap out old headers to match indexes where headers have been replaced
        * large header handler: You already indexed, but your .sam headers are too large to convert to .bam? This tool will generate a stand-alone header file and mapping script you can use to fix that.  

## Author
Jamie Alumbaugh <br>
Centre for Palaeogenetics, Stockholm
 
## Ackowledgements
Special thanks to NBIS (National Bioinformatics Infrastructure Sweden) for their Snakemake BYOC (bring-your-own-code) Workshop of Spring 2026. 

The basic github framework of this pipeline is based on [this template](https://github.com/snakemake-workflows/snakemake-workflow-template).

This pipeline was developed with documentation search and coding assistance from [Snakemake Guru AI](https://snakemake.readthedocs.io/en/stable/#) and [Claude Sonnet 5](https://claude.ai/chat/).


## Introduction

**Streit-lab/enhancer_annotation_and_motif_analysis** is a bioinformatic analysis pipeline for identifying enhancers associated to genes of interest and screening for motif binding sites.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a portable, reproducible manner.

## Pipeline summary

1. Conditionally unzip genome (--fasta) and GTF (--gtf) files
2. Filter genes of interest (--gene_list) from GTF and filter gene biotype entries in GTF
3. Conditionally extend length of peaks (--peaks) by a given length (--extend_peaks)
4. Retrieve filtered peak fasta sequences
5. Calculate background base frequencies for motif screening
6. Identify motif binding sites in peaks ([`fimo`](https://meme-suite.org/meme/doc/fimo.html))
7. Annotate peak-motif file with nearby genes

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)).

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run Streit-lab/enhancer_annotation_and_motif_analysis -r main -profile test,docker --outdir output
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   - Typical command for Streit-lab/enhancer_annotation_and_motif_analysis analysis:

   ```bash
   nextflow run Streit-lab/enhancer_annotation_and_motif_analysis --fasta <FASTA_PATH_OR_URL> --gtf <GTF_PATH_OR_URL> --motif_matrix <MEME_MOTIF_FILE> --peaks_bed <PEAK_BED_FILE> -profile <docker/singularity/conda>
   ```

## Pipeline parameters

`--fasta`: *Required*. Path or URL to fasta file, can be gzipped.

`--gtf`: *Required*. Path or URL to GTF file, can be gzipped.

`--motif_matrix`: *Required*. Path to matrix file in [`meme`](https://meme-suite.org/meme/doc/meme-format.html) format. [`Example file`](https://github.com/Streit-lab/enhancer_annotation_and_motif_analysis/blob/main/test_data/six1_motifs.txt). 

`--peaks_bed`: *Required*. Path to peak file in BED format. Must contain four columns; chrom, start, end, peakid. [`Example file`](https://github.com/Streit-lab/enhancer_annotation_and_motif_analysis/blob/main/test_data/peaks.bed).

`--gene_ids`: *Required*. List of gene ids present in GTF to screen for enhancers and motifs. One gene id per line. [`Example file`](https://github.com/Streit-lab/enhancer_annotation_and_motif_analysis/blob/main/test_data/peaks.bed).

`--extend_peaks`: *Optional*. Number of bases by which to extend peaks (up and downstream). Default = 0.

`--enhancer_window`: *Optional*. Distance from TSS in GTF within which enhancers are screened. Default = 50000.

`--markov_background`: *Optional*. Markov background model used to define base frequencies for motif screening. This is calculated by default from the provided --fasta input.

`--fimo_pval`: *Optional*. p-value threshold used by FIMO for motif screening.

`--gtf_gene_name_col`: *Optional*. Entry in GTF corresponding to gene names. Default = 'gene_name'.

`--gtf_gene_id_col`: *Optional*. Entry in GTF corresponding to gene names. Default = 'gene_id'.

`--outdir`: *Optional*. Directory to output results to. Default = 'results'.



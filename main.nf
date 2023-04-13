#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/
if(params.debug) {log.info Headers.build_debug_param_summary(params, params.monochrome_logs)}

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include { GUNZIP as GUNZIP_FASTA }      from "$baseDir/modules/local/gunzip/main"
include { FILTER_GTF_GENE_LIST }        from "$baseDir/modules/local/filter_gtf_gene_list/main"
include { ANNOTATE_PEAKS_TO_GTF }      from "$baseDir/modules/local/annotate_peaks_to_gtf/main"
include { BEDTOOLS_GETFASTA }           from "$baseDir/modules/nf-core/bedtools/getfasta/main"
include { MEMESUITE_FIMO }              from "$baseDir/modules/local/memesuite_fimo/main"
include { ANNOTATE_MOTIF_HITS }         from "$baseDir/modules/local/annotate_motif_hits/main"

/*------------------------------------------------------------------------------------
Set channels
--------------------------------------------------------------------------------------*/

Channel
    .value(params.motif_matrix)
    .set{ch_motif_matrix}

Channel
    .value(params.peaks_bed)
    .set{ch_peak_bed}

Channel
    .value(params.gene_ids)
    .set{ch_gene_ids}

Channel
    .value(params.gtf)
    .set{ch_gtf}

/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow {

    /*
    * Uncompress genome fasta file if required
    */
    if (params.fasta.endsWith(".gz")) {
        ch_fasta    = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    // Filter gtf based on gene list
    FILTER_GTF_GENE_LIST( ch_gtf, ch_gene_ids )

    // Annotate peaks based on GTF (either closest gene or using sliding window)
    ANNOTATE_PEAKS_TO_GTF(ch_peak_bed, FILTER_GTF_GENE_LIST.out.gtf)

    // Get fasta sequences for peaks which are associated to genes in gene list
    BEDTOOLS_GETFASTA(ANNOTATE_PEAKS_TO_GTF.out.gtf, ch_fasta)

    MEMESUITE_FIMO(ch_motif_matrix, BEDTOOLS_GETFASTA.out.fasta)
    ANNOTATE_MOTIF_HITS(ANNOTATE_PEAKS_TO_GTF.out.tsv, MEMESUITE_FIMO.out)
}
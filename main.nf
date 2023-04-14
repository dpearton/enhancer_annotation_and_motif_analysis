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
include { GUNZIP as GUNZIP_GTF }        from "$baseDir/modules/local/gunzip/main"
include { EXTEND_PEAKS }                from "$baseDir/modules/local/extend_peaks/main"
include { FILTER_GTF_GENE_LIST }        from "$baseDir/modules/local/filter_gtf_gene_list/main"
include { ANNOTATE_PEAKS_TO_GTF }       from "$baseDir/modules/local/annotate_peaks_to_gtf/main"
include { BEDTOOLS_GETFASTA }           from "$baseDir/modules/nf-core/bedtools/getfasta/main"
include { MEMESUITE_FASTA_GET_MARKOV }  from "$baseDir/modules/local/memesuite_fasta_get_markov/main"
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

/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow {
    // Uncompress genome fasta file if required
    if (params.fasta.endsWith(".gz")) {
        ch_fasta    = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file( params.fasta )
    }

    // Uncompress genome gtf file if required
    if (params.fasta.endsWith(".gz")) {
        ch_gtf    = GUNZIP_GTF ( params.gtf ).gunzip
    } else {
        ch_gtf = file( params.gtf )
    }

    // Filter gtf based on gene list
    FILTER_GTF_GENE_LIST( ch_gtf, ch_gene_ids )

    // Conditionally extend peak length
    // Annotate peaks based on GTF (either closest gene or using sliding window)
    if (params.extend_peaks != 0) {
        ANNOTATE_PEAKS_TO_GTF( EXTEND_PEAKS( ch_peak_bed, params.extend_peaks ).bed, FILTER_GTF_GENE_LIST.out.gtf )
    } else {
        ANNOTATE_PEAKS_TO_GTF( ch_peak_bed, FILTER_GTF_GENE_LIST.out.gtf )
    }

    // Get fasta sequences for peaks which are associated to genes in gene list
    BEDTOOLS_GETFASTA( ANNOTATE_PEAKS_TO_GTF.out.bed, ch_fasta )

    if (!params.markov_background){
        ch_background = MEMESUITE_FASTA_GET_MARKOV( ch_fasta )
    } else {
        ch_background = params.markov_background
    }

    // Find motifs in peak sequences
    MEMESUITE_FIMO( ch_motif_matrix, BEDTOOLS_GETFASTA.out.fasta, ch_background )

    // Add gene annotations to motif hits
    ANNOTATE_MOTIF_HITS( ANNOTATE_PEAKS_TO_GTF.out.tsv, MEMESUITE_FIMO.out )
}
#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Validate input parameters in specialised library
WorkflowParams.initialise(params, log)

// Print all params
// WorkflowMain.initialise(workflow, params, log)

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/
if(params.debug) {log.info Headers.build_debug_param_summary(params, params.monochrome_logs)}

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include { GUNZIP as GUNZIP_FASTA }                          from "$baseDir/modules/local/gunzip/main"
include { GUNZIP as GUNZIP_GTF }                            from "$baseDir/modules/local/gunzip/main"
include { EXTEND_PEAKS }                                    from "$baseDir/modules/local/extend_peaks/main"
include { EXTRACT_GTF_TRANSCRIPTS }                         from "$baseDir/modules/local/extract_gtf_transcripts/main"
include { FILTER_GTF_GENE_LIST }                            from "$baseDir/modules/local/filter_gtf_gene_list/main"
include { ANNOTATE_PEAKS_TO_GTF }                           from "$baseDir/modules/local/annotate_peaks_to_gtf/main"
include { ANNOTATE_PEAKS_TO_GTF_CTCF }                      from "$baseDir/modules/local/annotate_peaks_to_gtf_ctcf/main"
include { EXTRACT_FLANKING_CTCF }                           from "$baseDir/modules/local/extract_flanking_ctcf/main"
include { BEDTOOLS_GETFASTA }                               from "$baseDir/modules/nf-core/bedtools/getfasta/main"
include { BEDTOOLS_SORT as BEDTOOLS_SORT_PEAKS }            from "$baseDir/modules/local/bedtools/sort/main.nf"
include { BEDTOOLS_SORT as BEDTOOLS_SORT_FLANKING_CTCF }    from "$baseDir/modules/local/bedtools/sort/main.nf"
include { SAMTOOLS_FAIDX }                                  from "$baseDir/modules/local/samtools/faidx/main.nf"
include { MEMESUITE_FASTA_GET_MARKOV }                      from "$baseDir/modules/local/memesuite_fasta_get_markov/main"
include { MEMESUITE_FIMO }                                  from "$baseDir/modules/local/memesuite_fimo/main"
include { ANNOTATE_MOTIF_HITS }                             from "$baseDir/modules/local/annotate_motif_hits/main"

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

    SAMTOOLS_FAIDX(ch_fasta)

    // Uncompress genome gtf file if required
    if (params.fasta.endsWith(".gz")) {
        ch_gtf    = GUNZIP_GTF ( params.gtf ).gunzip
    } else {
        ch_gtf = file( params.gtf )
    }
    
    EXTRACT_GTF_TRANSCRIPTS( ch_gtf )

    // Filter gtf based on gene list if it is provided, otherwise annotate peaks using all genes
    if (params.gene_ids){
        ch_filtered_gtf = FILTER_GTF_GENE_LIST( EXTRACT_GTF_TRANSCRIPTS.out.gtf, ch_gene_ids )
    } else {
        ch_filtered_gtf = EXTRACT_GTF_TRANSCRIPTS.out.gtf
    }

    // Conditionally extend peak length
    if (params.extend_peaks != 0) {
        ch_peaks_processed = EXTEND_PEAKS( ch_peak_bed, params.extend_peaks ).bed
    } else {
        ch_peaks_processed = ch_peak_bed
    }

    BEDTOOLS_SORT_PEAKS ( ch_peaks_processed )

    // Annotate peaks based on GTF (either closest gene or using sliding window) OR within flanking CTCF sites
    if (params.ctcf) {
        EXTRACT_FLANKING_CTCF ( BEDTOOLS_SORT_PEAKS.out.sorted, params.ctcf, SAMTOOLS_FAIDX.out.fai )
        BEDTOOLS_SORT_FLANKING_CTCF ( EXTRACT_FLANKING_CTCF.out.bed )
        ANNOTATE_PEAKS_TO_GTF_CTCF ( ch_peak_bed, FILTER_GTF_GENE_LIST.out.gtf, BEDTOOLS_SORT_FLANKING_CTCF.out.sorted )

        ch_peak_annotations_bed = ANNOTATE_PEAKS_TO_GTF_CTCF.out.bed
        ch_peak_annotations_tsv = ANNOTATE_PEAKS_TO_GTF_CTCF.out.tsv

    } else {
        ANNOTATE_PEAKS_TO_GTF( ch_peak_bed, FILTER_GTF_GENE_LIST.out.gtf )

        ch_peak_annotations_bed = ANNOTATE_PEAKS_TO_GTF.out.bed
        ch_peak_annotations_tsv = ANNOTATE_PEAKS_TO_GTF.out.tsv
    }

    
    if(params.run_motif_analysis) {
    // Get fasta sequences for peaks which are associated to genes in gene list
    BEDTOOLS_GETFASTA( ch_peak_annotations_bed, ch_fasta )

    if (!params.markov_background){
        ch_background = MEMESUITE_FASTA_GET_MARKOV( ch_fasta )
    } else {
        ch_background = params.markov_background
    }

    // Find motifs in peak sequences
    MEMESUITE_FIMO( ch_motif_matrix, BEDTOOLS_GETFASTA.out.fasta, ch_background )

    // Add gene annotations to motif hits
    ANNOTATE_MOTIF_HITS( ch_peak_annotations_tsv, MEMESUITE_FIMO.out )
    }
}
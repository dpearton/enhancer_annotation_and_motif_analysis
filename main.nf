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
WorkflowMain.initialise(workflow, params, log)

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include { GUNZIP as GUNZIP_FASTA }                          from "$baseDir/modules/local/gunzip/main"
include { GUNZIP as GUNZIP_GTF }                            from "$baseDir/modules/local/gunzip/main"
include { GUNZIP as GUNZIP_GFF }                            from "$baseDir/modules/local/gunzip/main"
include { GFFREAD }                                         from "$baseDir/modules/nf-core/gffread/main"
include { EXTEND_PEAKS }                                    from "$baseDir/modules/local/extend_peaks/main"
include { EXTRACT_GTF_TRANSCRIPTS }                         from "$baseDir/modules/local/extract_gtf_transcripts/main"
include { FILTER_GTF_GENE_LIST }                            from "$baseDir/modules/local/filter_gtf_gene_list/main"
include { EXTRACT_GTF_WINDOW_COORDINATES }                  from "$baseDir/modules/local/extract_gtf_window_coordinates/main"
include { INTERSECT_PEAKS_GTF }                             from "$baseDir/modules/local/intersect_peaks_gtf/main"
include { COLLAPSE_BEDTOOLS_INTERSECT }                     from "$baseDir/modules/local/collapse_bedtools_intersect/main"
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

// Set motif_matrix to jaspar non redundant database if not provided
motif_matrix = params.motif_matrix ? params.motif_matrix : 'jaspar_core_vert_nonredundant_motifs'

if(motif_matrix == 'jaspar_core_vert_nonredundant_motifs') {
    ch_motif_matrix = Channel.from( file("$projectDir/assets/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt", checkIfExists: true))
} 
else if(motif_matrix == 'jaspar_core_vert_redundant_motifs') {
    ch_motif_matrix = Channel.from( file("$projectDir/assets/JASPAR2022_CORE_vertebrates_redundant_pfms_meme.txt", checkIfExists: true))
} 
else {
    ch_motif_matrix = Channel.from( file(motif_matrix, checkIfExists:true ))
}

Channel
    .from(file(params.peaks_bed, checkIfExists:true))
    .set{ch_peak_bed}

if (params.gene_ids){
Channel
    .from(file(params.gene_ids, checkIfExists:true))
    .set{ch_gene_ids}
}
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

    // Uncompress GTF annotation file or create from GFF3 if required
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( params.gtf ).gunzip
        } else {
            ch_gtf      = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff      = file(params.gff)
        }
        ch_gtf          = GFFREAD ( ch_gff ).gtf
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
        ANNOTATE_PEAKS_TO_GTF_CTCF ( ch_peak_bed, ch_filtered_gtf, BEDTOOLS_SORT_FLANKING_CTCF.out.sorted )

        ch_peak_annotations_bed = ANNOTATE_PEAKS_TO_GTF_CTCF.out.bed
        ch_peak_annotations_tsv = ANNOTATE_PEAKS_TO_GTF_CTCF.out.tsv

    } else {
        EXTRACT_GTF_WINDOW_COORDINATES( ch_filtered_gtf, SAMTOOLS_FAIDX.out.fai )
        INTERSECT_PEAKS_GTF( ch_peak_bed, EXTRACT_GTF_WINDOW_COORDINATES.out.bed )
        COLLAPSE_BEDTOOLS_INTERSECT( INTERSECT_PEAKS_GTF.out.bed )

        // ANNOTATE_PEAKS_TO_GTF( ch_peak_bed, ch_filtered_gtf )

        ch_peak_annotations_bed = COLLAPSE_BEDTOOLS_INTERSECT.out.bed
        ch_peak_annotations_tsv = COLLAPSE_BEDTOOLS_INTERSECT.out.tsv
    }

    
    if(!params.skip_motif_analysis) {
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
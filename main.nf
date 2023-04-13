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
include {MEMESUITE_FIMO} from "$baseDir/modules/local/memesuite_fimo/main"
include {BEDTOOLS_GETFASTA} from "$baseDir/modules/nf-core/bedtools/getfasta/main"
include {ANNOTATE_MOTIF_HITS} from "$baseDir/modules/local/annotate_motif_hits/main"

/*------------------------------------------------------------------------------------
Set channels
--------------------------------------------------------------------------------------*/

// Set channel for binary knowledge matrix for cell state classification
Channel
    .value(params.fasta)
    .set{ch_fasta}

Channel
    .value(params.motif_matrix)
    .set{ch_motif_matrix}

Channel
    .value(params.enhancer_bed)
    .set{ch_enhancer_bed}

Channel
    .value(params.annotated_enhancers)
    .set{ch_annotated_enhancers}

/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow {
    BEDTOOLS_GETFASTA(ch_enhancer_bed, ch_fasta)
    MEMESUITE_FIMO(ch_motif_matrix, BEDTOOLS_GETFASTA.out.fasta)
    ANNOTATE_MOTIF_HITS(ch_annotated_enhancers, MEMESUITE_FIMO.out)
}
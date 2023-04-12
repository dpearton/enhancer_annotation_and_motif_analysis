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

/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow {
    MEMESUITE_FIMO(ch_motif_matrix, ch_fasta)
}
//
// This file holds several functions specific to the Streit-lab/enhancer_annotation_and_motif_analysis pipeline
//

class WorkflowParams {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.gtf && !params.gff) {
            log.error "No GTF or GFF3 annotation specified! The pipeline requires at least one of these files."
            System.exit(1)
        }

        if (params.gtf && params.gff) {
            gtfGffWarn(log)
        }

        if (!params.peaks_bed) {
            log.error "Peaks bed file not specified with e.g. '--peaks_bed' or via a detectable config file."
            System.exit(1)
        }

        if (!params.gene_ids) {
            geneIdsWarn(log)
        }

        if (params.run_motif_analysis && !params.motif_matrix) {
            motifWarn(log)
        }
    }

    //
    // Print a warning if both GTF and GFF have been provided
    //
    private static void gtfGffWarn(log) {
        log.warn "=============================================================================\n" +
            "  Both '--gtf' and '--gff' parameters have been provided.\n" +
            "  Using GTF file as priority.\n" +
            "==================================================================================="
    }

    //
    // Print a warning if macs_gsize parameter has not been provided
    //
    private static void geneIdsWarn(log) {
        log.warn "=============================================================================\n" +
            "  --gene_ids parameter has not been provided.\n" +
            "  All gene IDs will be extracted from the genome annotation file and will be used to filter annotated peaks.\n" +
            "  This may result in the pipeline being extremely slow to run!\n" +
            "==================================================================================="
    }

    //
    // Print a warning if macs_gsize parameter has not been provided
    //
    private static void motifWarn(log) {
        log.warn "=============================================================================\n" +
            "  --motif_matrix parameter has not been provided.\n" +
            "  By default, the JASPAR core vertebrate redundant pfms collection will be used.\n" +
            "==================================================================================="
    }
}

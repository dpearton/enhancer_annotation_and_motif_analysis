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

        if (!params.gtf) {
            log.error "No GTF annotation specified!"
            System.exit(1)
        }

        if (!params.peaks_bed) {
            log.error "Peaks bed file not specified with e.g. '--peaks_bed' or via a detectable config file."
            System.exit(1)
        }

        if (params.run_motif_analysis) {
            if(!params.motif_matrix){
                log.error "No motif matrix provided with e.g. '--motif_matrix motifs.txt' or via a detectable config file."
                System.exit(1)
            }
        }
    }
}

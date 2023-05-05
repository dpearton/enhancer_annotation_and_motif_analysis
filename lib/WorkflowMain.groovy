//
// This file holds several functions specific to the main.nf workflow in the nf-core/cutandrun pipeline
//

class WorkflowMain {

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --fasta <FASTA_PATH_OR_URL> --gtf <GTF_PATH_OR_URL> --motif_matrix <MEME_MOTIF_FILE> --peaks_bed <PEAK_BED_FILE> -profile <docker/singularity/conda>"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }
}

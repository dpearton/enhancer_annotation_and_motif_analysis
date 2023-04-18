process ANNOTATE_PEAKS_TO_GTF_CTCF {
    label 'process_low'

    container "alexthiery/enhancer_annotation_and_motif_analysis:latest"

    input:
    path peaks
    path gtf
    path flanking_ctcf

    output:
    path "annotated_peaks.bed", emit: bed
    path "annotated_peaks.tsv", emit: tsv

    script:

    script:
    def args = task.ext.args  ?: ''

    """
    Rscript $baseDir/bin/annotate_peaks_to_gtf.R \\
        --peak_bed ${peaks} \\
        --gtf ${gtf} \\
        $args
    """
}
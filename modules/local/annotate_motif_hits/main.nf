process ANNOTATE_MOTIF_HITS {
    label 'process_low'

    container "alexthiery/enhancer_annotation_and_motif_analysis:latest"

    input:
    path annotated_peaks
    path motif_tsv

    output:
    path "annotated_motifs.tsv"

    script:
    def args = task.ext.args  ?: ''

    """
    Rscript $baseDir/bin/annotate_motif_hits.R \\
        --annotated_peaks ${annotated_peaks} \\
        --motif_tsv ${motif_tsv} \\
        $args
    """
}
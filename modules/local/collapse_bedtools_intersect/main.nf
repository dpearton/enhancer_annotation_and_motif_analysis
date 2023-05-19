process COLLAPSE_BEDTOOLS_INTERSECT {
    label 'process_low'

    container "alexthiery/enhancer_annotation_and_motif_analysis:latest"

    input:
    path peaks

    output:
    path "annotated_peaks.bed", emit: bed
    path "annotated_peaks.tsv", emit: tsv

    script:

    script:

    """
    Rscript $baseDir/bin/collapse_bedtools_intersect.R \\
        --peak_bed ${peaks}
    """
}
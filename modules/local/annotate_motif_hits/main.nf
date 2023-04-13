process ANNOTATE_MOTIF_HITS {
    label 'process_low'

    container "alexthiery/tidyverse_getopt:latest"

    input:
    path annotated_enhancers
    path motif_hits

    output:
    path "annotated_motif_hits.txt"

    script:
    def args = task.ext.args  ?: ''

    """
    Rscript $moduleDir/bin/annotate_motif_hits.R \\
        $args
    """
}
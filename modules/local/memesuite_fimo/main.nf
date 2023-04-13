process MEMESUITE_FIMO {
    label 'process_low'

    container "quay.io/biocontainers/meme:5.4.1--py310pl5321hb021246_2"

    input:
    path motif_matrix
    path fasta

    output:
    path "*"

    script:
    def args = task.ext.args  ?: ''


    """
    fimo \\
        $args \\
        ${motif_matrix} \\
        ${fasta}
    """
}
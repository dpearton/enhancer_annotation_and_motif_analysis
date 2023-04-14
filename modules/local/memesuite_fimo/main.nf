process MEMESUITE_FIMO {
    label 'process_low'

    container "quay.io/biocontainers/meme:5.4.1--py310pl5321hb021246_2"

    input:
    path motif_matrix
    path fasta
    path markov_background

    output:
    path "*"

    script:
    def args = task.ext.args  ?: ''


    """
    fimo \\
        $args \\
        --bfile ${markov_background} \\
        ${motif_matrix} \\
        ${fasta}
    """
}
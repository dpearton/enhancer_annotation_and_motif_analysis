process MEMESUITE_FASTA_GET_MARKOV {
    label 'process_low'

    container "quay.io/biocontainers/meme:5.4.1--py310pl5321hb021246_2"

    input:
    path fasta

    output:
    path "markov_background.txt"

    script:
    def args = task.ext.args  ?: ''

    """
    fasta-get-markov \\
        $args \\
        ${fasta} \\
        markov_background.txt
    """
}
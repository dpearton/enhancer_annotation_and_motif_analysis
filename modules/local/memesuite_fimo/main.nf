process MEMESUITE_FIMO {
    label 'process_low'

    container "memesuite/memesuite"

    input:
    path motif_matrix
    path fasta

    output:
    path "*.loom", emit: loom

    script:
    def args = task.ext.args  ?: ''


    """
    fimo \\
        $args \\
        ${motif_matrix} \\
        ${fasta}
    """
}
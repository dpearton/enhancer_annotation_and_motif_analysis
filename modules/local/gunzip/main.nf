process GUNZIP {
    tag "$archive"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path archive

    output:
    path "$gunzip",       emit: gunzip

    script:
    def args = task.ext.args ?: ''
    gunzip = archive.toString() - '.gz'
    """
    gunzip \\
        -f \\
        $args \\
        $archive
    """
}

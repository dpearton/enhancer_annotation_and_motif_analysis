process EXTRACT_GTF_WINDOW_COORDINATES {

    label 'process_high'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path(gtf_file)
    path(index_file)

    output:
    path "$output"       , emit: bed

    script:
    def args = task.ext.args  ?: ''

    output = gtf_file.toString() - ".gtf" + ".bed"
    
    """
    extract_gtf_window_coordinates.sh -g $gtf_file -i $index_file -o $output $args
    """
}
process INTERSECT_PEAKS_GTF {

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"
    
    input:
    path(peaks)
    path(gtf)

    output:
    path "$output"       , emit: bed

    script:
    output = peaks.toString() - ".bed" + "_peaks_intersected.bed"
    """
    bedtools intersect -a $peaks -b $gtf -wa -wb > $output
    """
}
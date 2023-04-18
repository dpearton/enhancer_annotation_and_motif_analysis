process EXTRACT_FLANKING_CTCF {
    label 'process_medium'

    container "alexthiery/parallel:latest"

    input:
    path peaks_bed
    path ctcf_bed
    path fasta_fai

    output:
    path "*.bed", emit: bed

    script:
    def prefix =  task.ext.prefix ?: 'ctcf_flanking'

    """
    parallel_ctcf_flanking_peaks.sh $peaks_bed $ctcf_bed $fasta_fai ${prefix}_peaks.bed
    """
}
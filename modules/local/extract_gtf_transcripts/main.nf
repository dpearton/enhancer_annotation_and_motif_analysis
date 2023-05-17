process EXTRACT_GTF_TRANSCRIPTS {
    label 'process_low'

    container "nfcore/base"

    input:
    path gtf

    output:
    path "*_transcript_extract.gtf", emit: gtf

    script:
    def prefix =  task.ext.prefix ?: gtf.toString() - '.gtf'

    """
    awk -F'\t' '\$3 ~ /transcript|mRNA/' ${gtf} > ${prefix}_transcript_extract.gtf
    """
}
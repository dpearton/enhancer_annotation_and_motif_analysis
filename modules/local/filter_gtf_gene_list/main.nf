process FILTER_GTF_GENE_LIST {
    label 'process_low'

    container "nfcore/base"

    input:
    path gtf
    path gene_ids

    output:
    path "*filtered.gtf", emit: gtf

    script:
    def prefix =  task.ext.prefix ?: gtf.toString() - '.gtf'

    """
    awk -F'\t' '\$3 ~ /CDS/' ${gtf} > temp.gtf
    grep -f ${gene_ids} temp.gtf > ${prefix}_filtered.gtf
    """
}
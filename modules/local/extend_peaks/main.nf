process EXTEND_PEAKS {
    label 'process_low'

    container "nfcore/base"

    input:
    path peaks
    val window

    output:
    path "*.bed", emit: bed

    script:
    def prefix =  task.ext.prefix ?: 'extended_'

    """
    awk -F'\t' -v OFS='\t' '{ \$2 = (\$2 - ${window} < 0 ? 0 : \$2 - ${window}); \$3 = \$3 + ${window}; print }' ${peaks} > ${prefix}${peaks}
    """
}
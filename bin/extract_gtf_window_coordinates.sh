#!/bin/bash

# Default values
window=50000
gene_id_key="gene_id"
gene_name_key="gene_name"

# Usage function
usage() {
    echo "Usage: $(basename "$0") -g <gtf_file> -i <index_file> -o <output_file> [-w <window>] [-x <gene_id_key>] [-n <gene_name_key>]"
    echo "Options:"
    echo "  -g <gtf_file>         GTF file"
    echo "  -i <index_file>       Index file"
    echo "  -o <output_file>      Output file"
    echo "  -w <window>           Window size (default: 50000)"
    echo "  -x <gene_id_key>     Gene ID key (default: gene_id)"
    echo "  -n <gene_name_key>   Gene Name key (default: gene_name)"
    exit 1
}

# Parse command line options
while getopts "g:i:o:w:x:n:" opt; do
    case $opt in
    g)
        gtf_file=$OPTARG
        ;;
    i)
        index_file=$OPTARG
        ;;
    o)
        output_file=$OPTARG
        ;;
    w)
        window=$OPTARG
        ;;
    x)
        gene_id_key=$OPTARG
        ;;
    n)
        gene_name_key=$OPTARG
        ;;
    \?)
        usage
        ;;
    esac
done

# Check required options
if [[ -z $gtf_file || -z $index_file || -z $output_file ]]; then
    usage
fi

while IFS=$'\t' read -r -a line; do
    # Extract the chromosome and store it as a variable
    chrom="${line[0]}"

    # Extend window from the gene start site depending on strand
    if [ "${line[6]}" == "+" ]; then
        start=$(( ${line[3]} - window ))
        end=$(( ${line[3]} + window ))
    elif [ "${line[6]}" == "-" ]; then
        start=$(( ${line[4]} - window ))
        end=$(( ${line[4]} + window ))
    fi

    # Extract the chromosome length from the index file
    chrom_length=$(awk -v chrom="$chrom" '$1 == chrom {print $2}' "$index_file")

    # Use the chromosome length to adjust the promoter coordinates
    if (( start < 0 )); then
        start=0
    fi
    if (( end > chrom_length )); then
        end=$chrom_length
    fi

    gene_id=$(echo "${line[8]}" | grep -oP "${gene_id_key} \"\K[^\"]+")
    gene_name=$(echo "${line[8]}" | grep -oP "${gene_name_key} \"\K[^\"]+")

    if [ -z "$gene_name" ]; then
        gene_name="$gene_id"
    fi

    # Print out the elements of interest
    echo -e "$chrom\t$start\t$end\t${line[5]}\t${line[6]}\t$gene_id\t$gene_name"

done < "$gtf_file" > temp.bed


# Remove duplicate lines from a file (without adjacent duplicates)
sort temp.bed | uniq > "$output_file"

rm temp.bed
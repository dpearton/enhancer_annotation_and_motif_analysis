#!/bin/bash

# Define input file names
peaks_file="$1"
ctcf_file="$2"
fa_fai_file="$3"
output_file="$4"

# Calculate the number of lines in the input file
total_lines=$(wc -l < "$peaks_file")
lines_per_chunk=$((total_lines/10))

# Split the input file into chunks
split -l "$lines_per_chunk" "$peaks_file" "$peaks_file."

# Process each chunk in parallel
for file in "$peaks_file."*; do
  {
    # Loop through each row in the chunk and write output to a new bed file
    while read line; do

      # Extract the peak chromosome and store as a variable
      peak_chrom=$(echo "$line" | awk '{print $1}')

      # Extract the peak start and store as a variable
      peak_start=$(echo "$line" | awk '{print $2}')

      # Extract the peak end and store as a variable
      peak_end=$(echo "$line" | awk '{print $3}')

      # Extract the peak end and store as a variable
      peak_id=$(echo "$line" | awk '{print $4}')

      # Filter ctcf.bed by chromosome and capture output in a variable
      peak_match=$(awk -F '\t' -v chrom="$peak_chrom" '$1 == chrom' "$ctcf_file")

      # Store the closest but smaller value in a variable
      closest_smaller=$(echo "$peak_match" | awk -v start="$peak_start" '$3 < start {print $3}' | tail -n 1)

      # Store the closest but greater value in a variable
      closest_greater=$(echo "$peak_match" | awk -v end="$peak_end" '$2 > end {print $2}' | head -n 1)

      # Some peaks do not have flanking CTCF sites - i.e. if they are at the start or the end of the chromosome
      # If closest_smaller is empty, then set the value as 0 as this is the start of the chromosome
      if [[ -z "$closest_smaller" ]]; then
        closest_greater=0
      fi

      # If closest_greater is empty, use the peak_chrom variable to select the row from the fa.fai file
      # and extract the chromosome sequence length from the second column
      if [[ -z "$closest_greater" ]]; then
        chrom_length=$(awk -v chrom="$peak_chrom" '$1 == chrom {print $2}' "$fa_fai_file")
        closest_greater="$chrom_length"
      fi

      # Write the output to a temporary file
      echo -e "${peak_chrom}\t${closest_smaller}\t${closest_greater}\t${peak_id}" >> "output.$$.$file.tmp"

    done < "$file"

    # Concatenate all temporary files into a single output file
    cat "output.$$.$file.tmp" >> "output.$$"

    # Remove the temporary file
    rm "output.$$.$file.tmp"

  } &
done

# Wait for all processes to finish
wait

# Move the concatenated output file to the final output file
mv "output.$$" "$output_file"

# Remove the split input files
rm "$peaks_file."*

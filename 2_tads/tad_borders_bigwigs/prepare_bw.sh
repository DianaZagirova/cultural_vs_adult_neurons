#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

input_dir="$1"

# Path to the chrom.sizes file
chrom_sizes_path="../4_additional_files/chrom.sizes"

if [ ! -f "$chrom_sizes_path" ]; then
  echo "Error: Chrom.sizes file not found at $chrom_sizes_path"
  exit 1
fi

for bedgraph_file in "$input_dir"/*.bed; do
  # Check if the file exists (in case no .bedGraph files are found)
  if [ ! -e "$bedgraph_file" ]; then
    echo "No .bed files found in the directory."
    exit 1
  fi

  output_file="${bedgraph_file%.bed}.bw"

  bedGraphToBigWig "$bedgraph_file" "$chrom_sizes_path" "$output_file"

  if [ $? -eq 0 ]; then
    echo "Successfully converted $bedgraph_file to $output_file"
  else
    echo "Error converting $bedgraph_file"
  fi
done
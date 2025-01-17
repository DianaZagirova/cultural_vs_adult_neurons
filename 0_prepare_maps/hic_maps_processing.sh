#!/bin/bash

# Master script for HiC data processing: merge samples, drop diagonal, downsamaple, zoomify and balance
# Usage: ./master_hic_processing.sh -o OUTPUT_DIR -d DOWNSAMPLE_VALUE -g GENOME -f FILE1 [FILE2 ...]
# Example commands
## For .cool files
##./master_hic_processing.sh -o output_dir -d 1000000 -f file1.cool file2.cool

## For .mcool files (will automatically use 5kb resolution)
##./master_hic_processing.sh -o output_dir -d 1000000 -f file1.mcool file2.mcool

# Default values
OUTPUT_DIR=""
DOWNSAMPLE_VALUE=""
TEMP_DIR=""
GENOME="hg38"  # Default genome assembly
declare -a INPUT_FILES

# Parse command line arguments
while getopts "o:d:g:f:" opt; do
  case $opt in
    o) OUTPUT_DIR="$OPTARG"
       ;;
    d) DOWNSAMPLE_VALUE="$OPTARG"
       ;;
    g) GENOME="$OPTARG"
       ;;
    f) # Handle multiple input files
       while [[ "$#" -gt 0 ]]; do
         case $1 in
           -o|-d|-g|-f) break ;;
           *) INPUT_FILES+=("$1") ;;
         esac
         shift
       done
       ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
        ;;
  esac
done

# Validate required parameters
if [ -z "$OUTPUT_DIR" ] || [ -z "$DOWNSAMPLE_VALUE" ] || [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "Usage: $0 -o OUTPUT_DIR -d DOWNSAMPLE_VALUE [-g GENOME] -f FILE1 [FILE2 ...]"
    echo "  -o: Output directory"
    echo "  -d: Downsample value (integer)"
    echo "  -g: Genome assembly (default: hg38, also supports mm10, mm39)"
    echo "  -f: Input cool files (one or more)"
    exit 1
fi

# Validate genome assembly
case $GENOME in
    "hg38"|"mm10"|"mm39") ;;
    *)  echo "Error: Unsupported genome assembly. Supported values: hg38, mm10, mm39"
        exit 1
        ;;
esac

# Create output directories
mkdir -p "$OUTPUT_DIR"/{temp,merged,no_diag,final}
TEMP_DIR="$OUTPUT_DIR/temp"
MERGED_DIR="$OUTPUT_DIR/merged"
NO_DIAG_DIR="$OUTPUT_DIR/no_diag"
FINAL_DIR="$OUTPUT_DIR/final"

echo "Starting HiC data processing pipeline..."

# Function to check if a cooler file is valid
check_cooler_file() {
    local file=$1
    local step=$2
    
    if [ ! -f "$file" ]; then
        echo "Error: File not created at $step step: $file"
        exit 1
    }
    
    # Check if the file is a valid cooler file
    if ! cooler info "$file" &>/dev/null; then
        echo "Error: Invalid cooler file created at $step step: $file"
        exit 1
    }
    
    echo "✓ Validation passed for $step: $(basename "$file")"
}

# Function to check file size
check_file_size() {
    local file=$1
    local min_size=1000  # 1KB minimum size
    
    local size=$(stat -f%z "$file")
    if [ "$size" -lt "$min_size" ]; then
        echo "Warning: File $file seems too small (${size} bytes)"
        exit 1
    }
}

# Step 1: Merge input files if more than one
MERGED_FILE="$MERGED_DIR/merged.cool"
if [ ${#INPUT_FILES[@]} -gt 1 ]; then
    echo "Merging ${#INPUT_FILES[@]} input files..."
    
    # Check if all input files are .mcool
    ALL_MCOOL=true
    for file in "${INPUT_FILES[@]}"; do
        if [[ ! $file =~ \.mcool$ ]]; then
            ALL_MCOOL=false
            break
        fi
        # Validate input files
        if [[ $file =~ \.mcool$ ]]; then
            if ! cooler info "${file}::/resolutions/5000" &>/dev/null; then
                echo "Error: Invalid or missing 5kb resolution in mcool file: $file"
                exit 1
            fi
        else
            check_cooler_file "$file" "input"
        fi
    done
    
    if [ "$ALL_MCOOL" = true ]; then
        echo "Detected .mcool files, merging at 5kb resolution..."
        MCOOL_FILES=()
        for file in "${INPUT_FILES[@]}"; do
            MCOOL_FILES+=("${file}::/resolutions/5000")
        done
        cooler merge "$MERGED_FILE" "${MCOOL_FILES[@]}"
        check_cooler_file "$MERGED_FILE" "merge"
    else
        echo "Merging .cool files..."
        cooler merge "$MERGED_FILE" "${INPUT_FILES[@]}"
        check_cooler_file "$MERGED_FILE" "merge"
    fi
    WORKING_FILE="$MERGED_FILE"
else
    WORKING_FILE="${INPUT_FILES[0]}"
    # If single input is mcool, extract 5kb resolution
    if [[ $WORKING_FILE =~ \.mcool$ ]]; then
        echo "Single .mcool file detected, extracting 5kb resolution..."
        cp "$MERGED_FILE" "${WORKING_FILE}::/resolutions/5000"
    else
        cp "$WORKING_FILE" "$MERGED_FILE"
    fi
fi

# Step 2: Drop diagonal
echo "Dropping diagonal..."
FILE_NAME=$(basename "$WORKING_FILE" .cool)

# Check if bins and pixels files are created and not empty
cooler dump -t bins "$WORKING_FILE" > "$TEMP_DIR/${FILE_NAME}.bins.txt"
check_file_size "$TEMP_DIR/${FILE_NAME}.bins.txt"

cooler dump -t pixels --chunksize 100000000 "$WORKING_FILE" | awk -F'\t' '$1!=$2' > "$TEMP_DIR/${FILE_NAME}.pixels.txt"
check_file_size "$TEMP_DIR/${FILE_NAME}.pixels.txt"

NO_DIAG_FILE="$NO_DIAG_DIR/${FILE_NAME}.no_diag.cool"
cooler load -f coo --assembly "$GENOME" --chunksize 100000000 \
    "$TEMP_DIR/${FILE_NAME}.bins.txt" \
    "$TEMP_DIR/${FILE_NAME}.pixels.txt" \
    "$NO_DIAG_FILE"

check_cooler_file "$NO_DIAG_FILE" "diagonal removal"

# Step 3: Downsample
echo "Downsampling to target value: $DOWNSAMPLE_VALUE..."
FINAL_FILE="$FINAL_DIR/${FILE_NAME}.final.cool"
cooltools random-sample -c "$DOWNSAMPLE_VALUE" --exact -p 19 --chunksize 100000000 \
    "$NO_DIAG_FILE" "$FINAL_FILE"

check_cooler_file "$FINAL_FILE" "downsampling"

# Create multi-resolution cooler file
echo "Creating multi-resolution file..."
MCOOL_FILE="${FINAL_FILE%.cool}.mcool"
cooler zoomify -p 18 --balance \
    -r 5000,10000,15000,20000,25000,50000,100000,150000,200000,250000,500000,1000000,2000000,5000000 \
    "$FINAL_FILE"

# Validate final mcool file
if ! cooler info "${MCOOL_FILE}::/resolutions/5000" &>/dev/null; then
    echo "Error: Failed to create valid multi-resolution file"
    exit 1
fi

echo "✓ All processing steps completed successfully!"
echo "Final output files:"
echo "- Cool file: $FINAL_FILE"
echo "- Multi-resolution file: $MCOOL_FILE"

# Print file sizes
echo -e "\nFile sizes:"
ls -lh "$FINAL_FILE" "$MCOOL_FILE" | awk '{print $9 ": " $5}'

# Cleanup temporary files
rm -f "$TEMP_DIR"/*.txt

echo "Processing complete!"


# Generalized Imaging File Finder Script
# 
# find ADNI-3/sub-ADNI021S6910/ses-1/anat -type f -name "*T1w.nii*"
#
# === USER DEFINED VARIABLES ===
# Debug flag: set to 1 for verbose output, 0 for minimal output
DEBUG=0
# Output directory (will be created if it doesn't exist)
OUTPUT_DIR="/space/mcdonald-syn01/1/data/ai/2026_06_01/00_reports/file_check/cmig_struct_func_session1_run1/"  # <-- Set this to your desired output directory

# List of datasets
DATASETS=( CAPES DLBS_dataset ECP_dataset_bids ENIGMA_dataset HEP_dataset_bids IDEAS_dataset MTL_proj VU_dataset ADNI-3)

# List of derivative sub-paths (relative to dataset root)
DERIV_PATHS=(
  "sub-*/ses*1*/anat+sub-*/ses-pre/anat+sub-*/anat"         # For T1w raw files
  "derivatives/cat12-r2166/sub-*/ses*1*+derivatives/cat12-r2166/sub-*/ses-pre+derivatives/cat12-r2166/sub-*" # For mwp1 processed files
  "sub-*/ses*1*/func+sub-*/ses-pre/func+sub-*/func" # For func files: accumulate ses-01, ses-pre, and no-session subjects
  "derivatives/fmri_postproc/sub-*/ses*1*/reho+derivatives/fmri_postproc/sub-*/ses-pre/reho+derivatives/fmri_postproc/sub-*/reho" # For postprocessed reho files
  "derivatives/fmri_postproc/sub-*/ses*1*/alff+derivatives/fmri_postproc/sub-*/ses-pre/alff+derivatives/fmri_postproc/sub-*/alff" # For postprocessed alff files
  "derivatives/fmri_postproc/sub-*/ses*1*/alff+derivatives/fmri_postproc/sub-*/ses-pre/alff+derivatives/fmri_postproc/sub-*/alff" # For postprocessed falff files
)

# List of file patterns (relative to derivative path)
FILE_PATTERNS=(
  "*T1w.nii*"         # For T1w raw files (.nii and .nii.gz)
  "anat/*mwp1sub*.nii*" # For mwp1 processed files (.nii and .nii.gz)
  "*run*1*_bold.nii*;*_bold.nii*"        # For func files (.nii and .nii.gz)
  "*run*1_*desc-postprocessed_reho_z_smooth-4mm.nii*;*desc-postprocessed_reho_z_smooth-4mm.nii*" # For postprocessed reho files (.nii and .nii.gz)
  "*run*1_*desc-postprocessed_ALFF_z_smooth-4mm.nii*;*desc-postprocessed_ALFF_z_smooth-4mm.nii*" # For postprocessed alff files (.nii and .nii.gz)
  "*run*1_*desc-postprocessed_fALFF_z_smooth-4mm.nii*;*desc-postprocessed_fALFF_z_smooth-4mm.nii*" # For postprocessed falff files (.nii and .nii.gz)
)

# List of output file name suffixes (will be prefixed with dataset name)
OUTPUT_SUFFIXES=(
  "_t1w_raw.txt"
  "_mwp1_processed.txt"
  "_func_raw.txt"
  "_reho_postproc.txt"
  "_alff_postproc.txt"
  "_falff_postproc.txt"
)

# Strings to exclude from results (any file path containing these substrings will be filtered out)
EXCLUDE_PATTERNS=("Hypercapnia" "Scenes" "VentralVisual" "Words" "scenes" "language" "semantic" "social" "PURE" )

# === SCRIPT LOGIC ===

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check that all arrays are the same length
if [ ${#DERIV_PATHS[@]} -ne ${#FILE_PATTERNS[@]} ] || [ ${#DERIV_PATHS[@]} -ne ${#OUTPUT_SUFFIXES[@]} ]; then
  echo "Error: DERIV_PATHS, FILE_PATTERNS, and OUTPUT_SUFFIXES arrays must have the same length." >&2
  exit 1
fi

# Main loop
for dataset in "${DATASETS[@]}"; do
  for i in "${!DERIV_PATHS[@]}"; do
    deriv_paths="${DERIV_PATHS[$i]}"
    file_patterns="${FILE_PATTERNS[$i]}"
    output_file="$OUTPUT_DIR/${dataset}${OUTPUT_SUFFIXES[$i]}"

    # Separator ';' = try paths in order, stop at first match (fallback)
    # Separator '+' = try all paths and accumulate results
    if [[ "$deriv_paths" == *"+"* ]]; then
      IFS='+' read -ra PATHS <<< "$deriv_paths"
      accumulate_all=1
    else
      IFS=';' read -ra PATHS <<< "$deriv_paths"
      accumulate_all=0
    fi
    IFS=';' read -ra PATTERNS <<< "$file_patterns"
    > "$output_file" # Empty output file before writing
    found_any=0
    for path in "${PATHS[@]}"; do
      search_path="$dataset/$path"
      if [ "$DEBUG" -eq 1 ]; then
        echo "Globbing and searching in: $search_path"
        echo "Patterns: ${PATTERNS[*]}"
        echo "Output: $output_file"
      fi
      found_in_this_path=0
      for expanded_path in $search_path; do
        if [[ -d "$expanded_path" ]]; then
          for pattern in "${PATTERNS[@]}"; do
            tmp_pat=$(mktemp)
            if [[ "$pattern" == *"/"* ]]; then
              # Use find with -path for patterns with slashes
              if [ "$DEBUG" -eq 1 ]; then
                echo "Running: find $expanded_path -type f -path '$expanded_path/$pattern'"
              fi
              find -L "$expanded_path" -type f -path "$expanded_path/$pattern" >> "$tmp_pat"
            else
              # Use shell globbing for patterns without slashes
              if [ "$DEBUG" -eq 1 ]; then
                echo "Running: for f in $expanded_path/$pattern; do [ -f \"$f\" ] && echo \"$f\"; done"
              fi
              for f in $expanded_path/$pattern; do
                [ -f "$f" ] && echo "$f" >> "$tmp_pat"
              done
            fi
            # Filter out any lines matching EXCLUDE_PATTERNS
            for excl in "${EXCLUDE_PATTERNS[@]}"; do
              grep -v "$excl" "$tmp_pat" > "${tmp_pat}.filtered" && mv "${tmp_pat}.filtered" "$tmp_pat"
            done
            pat_count=$(wc -l < "$tmp_pat")
            if [ "$pat_count" -gt 0 ]; then
              cat "$tmp_pat" >> "$output_file"
              found_in_this_path=1
              rm -f "$tmp_pat"
              break  # Stop after first pattern that yields results
            fi
            rm -f "$tmp_pat"
          done
        fi
      done
      count=$(wc -l < "$output_file")
      if [ "$count" -gt 0 ]; then
        found_any=1
        if [ "$accumulate_all" -eq 0 ]; then
          break  # ';' mode: stop at first path with results
        fi
      fi
    done
    echo "[$dataset] $output_file: $count files found."
    if [ "$count" -eq 0 ]; then
      echo "  Warning: No files found for $dataset/($deriv_paths)/($file_patterns)" >&2
    fi
  done
done






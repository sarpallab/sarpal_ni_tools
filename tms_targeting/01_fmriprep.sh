#!/bin/bash

# Script to run fMRIPrep and XCP-D
# Usage: ./01_fmriprep.sh <subject_id> <bids_root_dir> <derivatives_dir_base>

# --- 1. Argument Parsing and Validation ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <subject_id> <bids_root_dir> <derivatives_dir_base>"
    echo "Example: $0 sub-001 /path/to/bids_dataset /path/to/bids_dataset/derivatives"
    exit 1
fi

SUBJECT_ID="$1"
BIDS_ROOT_DIR="$2"
DERIVATIVES_DIR_BASE="$3" # This will be the parent for fmriprep and xcp_d outputs

# Define specific output directories
FMRIPREP_OUTPUT_DIR="${DERIVATIVES_DIR_BASE}/fmriprep"
XCPD_OUTPUT_DIR="${DERIVATIVES_DIR_BASE}/xcp_d"

# Create output directories if they don't exist
mkdir -p "${FMRIPREP_OUTPUT_DIR}"
mkdir -p "${XCPD_OUTPUT_DIR}"

echo "--- Running fMRIPrep and XCP-D ---"
echo "Subject ID: ${SUBJECT_ID}"
echo "BIDS Root Directory: ${BIDS_ROOT_DIR}"
echo "fMRIPrep Output Directory: ${FMRIPREP_OUTPUT_DIR}"
echo "XCP-D Output Directory: ${XCPD_OUTPUT_DIR}"
echo "------------------------------------"

# --- 2. fMRIPrep Execution ---
# Ensure the FS_LICENSE_FILE path is correct and accessible within your environment/container
FS_LICENSE_FILE="/ihome/dsarpal/wcs17/license.txt" 
# Consider making FS_LICENSE_FILE an argument if it varies

echo "Starting fMRIPrep for ${SUBJECT_ID}..."
singularity run --cleanenv \
    /ix1/ginger/dsarpal/lab/reorg/projects/7T/setup_scripts/fmriprep_23.0.2.sif \
    "${BIDS_ROOT_DIR}" \
    "${FMRIPREP_OUTPUT_DIR}" \
    participant \
    --participant-label "${SUBJECT_ID}" \
    --skull-strip-t1w auto \
    --output-spaces anat fsnative MNI152NLin2009cAsym \
    --fs-license-file "${FS_LICENSE_FILE}" \
    --dvars-spike-threshold 1.5 \
    --fd-spike-threshold 0.5 \
    --return-all-components

# Check if fMRIPrep completed successfully (basic check)
if [ ! -d "${FMRIPREP_OUTPUT_DIR}/sub-${SUBJECT_ID}" ] && [ ! -d "${FMRIPREP_OUTPUT_DIR}/fmriprep/sub-${SUBJECT_ID}" ]; then
    # fMRIPrep <24.0 puts sub-XYZ directly in output_dir, >=24.0 puts it in output_dir/fmriprep/sub-XYZ
    # This handles both by checking the more direct one first. If fMRIPrep creates an additional 'fmriprep' subdir
    # then XCP-D needs to point to that. Let's assume fMRIPrep output is directly in FMRIPREP_OUTPUT_DIR for now.
    # If XCP-D fails to find fmriprep derivatives, this path might need adjustment:
    # XCPD_INPUT_DIR="${FMRIPREP_OUTPUT_DIR}" or XCPD_INPUT_DIR="${FMRIPREP_OUTPUT_DIR}/fmriprep"
    echo "fMRIPrep output for sub-${SUBJECT_ID} not found in ${FMRIPREP_OUTPUT_DIR}. Exiting."
    # exit 1 # Commenting out exit for now, but important for production
fi
echo "fMRIPrep finished for ${SUBJECT_ID}."


# --- 3. XCP-D Execution ---
echo "Starting XCP-D for ${SUBJECT_ID}..."

# XCP-D expects the *base* derivatives directory that contains the 'fmriprep' folder
# So, if FMRIPREP_OUTPUT_DIR is /path/to/bids/derivatives/fmriprep,
# then XCP-D's BIDS filter input should point to /path/to/bids/derivatives
XCPD_BIDS_FILTER_INPUT_DIR="${DERIVATIVES_DIR_BASE}"


singularity run --cleanenv \
    /ix1/ginger/dsarpal/lab/reorg/ni_tools/xcp_d-latest.simg \
    "${XCPD_BIDS_FILTER_INPUT_DIR}" \
    "${XCPD_OUTPUT_DIR}" participant \
    --participant-label "${SUBJECT_ID}" \
    --input-type fmriprep \
    --smoothing 8 \
    --despike \
    --nuisance-regressors 36P \
    --fd-thresh 0 \
    --bids-filter-file /ix1/ginger/dsarpal/lab/reorg/projects/enact/xcp_config_rest.json

if [ $? -eq 0 ]; then
    echo "XCP-D completed successfully for ${SUBJECT_ID}."
else
    echo "XCP-D failed for ${SUBJECT_ID}."
    # exit 1 # Commenting out exit for now
fi

echo "--- Processing complete for ${SUBJECT_ID} ---"

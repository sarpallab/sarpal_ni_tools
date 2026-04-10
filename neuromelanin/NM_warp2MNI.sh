{\rtf1\ansi\ansicpg1252\cocoartf2869
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/bash\
# Defines the interpreter as bash.\
\
# ---------------------------------------------------------\
# 1. ENVIRONMENT & VARIABLE SETUP\
# ---------------------------------------------------------\
# Load required software modules (FSL for brain extraction, ANTs for registration, MATLAB)\
module load fsl ants matlab\
\
# Define base directory paths for the study and data\
Study_dir="/eru/shares/sarpal/data/SarpalLab/DARES/7TPLUS"\
data_dir="/eru/shares/sarpal/data/SarpalLab/DARES/7TPLUS"\
\
# (Commented out) Control data directory\
#ctrl_dir="/ix1/ginger/dsarpal/lab/reorg/projects/neurolemanin/Z10_control_nifti"\
\
# Define paths to standard templates\
# Standard MNI152 template (1mm resolution) used as the target for registration\
MNI_template="/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/scripts/MNI152_T1_1mm_brain.nii.gz"\
# Study-specific Neuromelanin template\
NM_template="/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/Pitt_NM/Pitt_NM_template0.nii.gz"\
\
# (Commented out) Block to check for and create an "all_NM" directory if it doesn't exist\
#if; then\
#   echo "all_NM directory already exist"\
#else \
#   mkdir -p "$\{Study_dir\}/all_NM"\
#fi\
\
\
# ---------------------------------------------------------\
# 2. SUBJECT & SESSION LOOPS\
# ---------------------------------------------------------\
# NOTE: Currently, this loop is hardcoded to target exactly ONE subject (sub-2927). \
# If you want to loop through all subjects, you should use: for sub in "$\{data_dir\}"/sub-*; do\
for sub in "$\{data_dir\}"/sub-2927; do\
    # Extract just the folder name (e.g., "sub-2927") from the full path\
    subj=$(basename "$sub")\
\
    # NOTE: Hardcoded to exactly ONE session (ses-20251212).\
    # To loop all sessions, use: for ses in "$\{sub\}"/ses-*; do\
    for ses in "$\{sub\}"/ses-20251212; do\
        # Extract just the session name\
        sesh=$(basename "$ses")\
\
        # ---------------------------------------------------------\
        # 3. BRAIN EXTRACTION (BET)\
        # ---------------------------------------------------------\
        # Check if the skull-stripped brain file already exists\
        if; then\
            # If it exists, skip to save time. ($0 prints the name of the script)\
            echo "$\{0\}: $\{subj\}_T1w_BET.nii.gz file already exists"\
        else\
            # If it doesn't exist, run FSL's bet2 (Brain Extraction Tool)\
            # Input: T1w_INV2 image\
            # Output: T1w_BET image\
            # -f 0.4: Fractional intensity threshold (smaller values give larger brain outlines)\
            bet2 $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_T1w_INV2.nii.gz $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_T1w_BET.nii.gz -f 0.4\
        fi\
\
        # ---------------------------------------------------------\
        # 4. REGISTRATION TO MNI SPACE\
        # ---------------------------------------------------------\
        # Check if the registered T1 file (warped) already exists\
        if; then\
            echo "$\{0\}: $\{subj\}_rT1w.nii.gz file already exists"\
        else\
            # Register the skull-stripped T1w image to the MNI template using ANTs\
            # -d 3: 3-Dimensional image\
            # -f: Fixed image (the MNI template)\
            # -m: Moving image (the subject's T1w brain)\
            # -o: Output prefix (generates Warped, InverseWarped, and transform matrices)\
            antsRegistrationSyNQuick.sh -d 3 \\\
                -f $\{MNI_template\} \\\
                -m $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_T1w_BET.nii.gz \\\
                -o $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_rT1w_BET_\
        fi\
\
        # ---------------------------------------------------------\
        # 5. NEUROMELANIN BIAS CORRECTION & REGISTRATION\
        # ---------------------------------------------------------\
        # BUG WARNING: You are checking if a file exists in the "nm_data" directory, \
        # but your commands below save the output to the "anat" directory. \
        # You may want to fix this path discrepancy.\
        if; then\
            echo "$\{0\}: $\{subj\}_run-1_rNM.nii.gz file already exists"\
        else\
            # Apply N4 Bias Field Correction to the Neuromelanin image to fix intensity unevenness\
            N4BiasFieldCorrection -d 3 \\\
                -i $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_acq-NMMT_MTR.nii.gz \\\
                -o $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_bias_nm.nii.gz\
            \
            # Apply the transform matrices calculated during T1->MNI registration \
            # to move the bias-corrected NM image into MNI space\
            # -r: Reference image (MNI space)\
            # -t: Transform files (applied in order: Warp field first, then Affine matrix)\
            antsApplyTransforms -d 3 \\\
                -i $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_bias_nm.nii.gz \\\
                -r $\{MNI_template\} \\\
                -o $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_bias_rnm.nii.gz \\\
                -t $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_rT1w_BET_1Warp.nii.gz \\\
                -t $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_rT1w_BET_0GenericAffine.mat\
            \
            # Unzip the final registered Neuromelanin NIfTI file (removes the .gz extension)\
            gunzip $\{data_dir\}/$\{subj\}/$\{sesh\}/anat/$\{subj\}_$\{sesh\}_bias_rnm.nii.gz \
        fi\
\
    done\
done}
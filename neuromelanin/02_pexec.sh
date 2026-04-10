#! /bin/bash
# Defines the interpreter as bash.

# ---------------------------------------------------------
# 1. SETUP DIRECTORIES & INPUT FILES
# ---------------------------------------------------------
# Define the directory containing the processed NM images
study_dir="/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/Pitt_NM/NM"

# Define the directory where the temporary and final template files will be saved
out_dir="/eru/shares/sarpal/data/SarpalLab/DARES/7T/Neuromelanin/Pitt_NM/out_temp"

# Create a bash array named 'images' containing all unzipped .nii files in the study directory
images=($study_dir/*.nii)

# ---------------------------------------------------------
# 2. RUN TEMPLATE CONSTRUCTION
# ---------------------------------------------------------
# Execute the ANTs script to build an unbiased population template
# -d 3: Specifies 3-Dimensional images
# -o $out_dir/Pitt_NM_: Specifies the output directory and the prefix for all generated files
# -c 2: Runs the SyN (Symmetric Normalization) registration algorithm (Type 2 in ANTs)
# -j 8: Utilizes 8 CPU cores for parallel processing to speed up the registration
# "${images}": Expands the array of images to pass every .nii file as an input argument
./antsMultivariateTemplateConstruction2.sh \
    -d 3 \
    -o $out_dir/Pitt_NM_ \
    -c 2 \
    -j 8 \
    "${images}"
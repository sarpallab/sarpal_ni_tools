template="/ihome/crc/install/fsl/6.0.4/centos7/fsl/data/standard/MNI152_T1_1mm.nii.gz"
#left_mask="left_hemisphere_mask.nii.gz"
#right_mask="right_hemisphere_mask.nii.gz"

#for folder in /ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/sub-2346/ses-20181221/; do
for folder in /ix1/ginger/dsarpal/lab/reorg/projects/mrsi/Z5_mrsi/SCZ/sub*/ses*/; do

	ref=$(basename $template)

        ln -sf ${template} ${folder}/${ref}
	
	sub=$(basename $(dirname $folder))
	ses=$(basename $folder)
	prefix=${sub}_${ses}

	scout_prefix=${prefix}"_scout_tonative_"
	native_prefix=${prefix}"_native_tomni_"


	TARGET_REGIONS=("lCaud" "rCaud" "lThal" "rThal")

	for region in "${TARGET_REGIONS[@]}"; do
    
		echo
		echo "Processing region: ${region} for ${prefix}"

		# --- 1. Prepare region-specific input point file ---
		# This file will contain original x,y,z,t from the master file for the current SUBJ, SES, and ROI.
		MASTER_POINTS_FILE_PATH="./points.csv"
		# Columns in this file will be: x,y,z,t (comma-separated)
		region_specific_input_points_file="${folder}/points_input_for_${region}.csv"
		rm $region_specific_input_points_file
    
		# Extract x,y,z,t (columns 3,4,5,6 from master) for the current SUBJ, SES, and ROI label (column 7 from master)
		awk_comparison_subj="${sub#sub-}"
		awk_comparison_ses="${ses#ses-}"

	tr -d '\r' < "${MASTER_POINTS_FILE_PATH}" | awk -F, -v subj="${awk_comparison_subj}" -v ses="${awk_comparison_ses}" -v label="${region}" \
		'NR > 1 && $1 == subj && $2 == ses && $7 == label {OFS=","; print $3,$4,$5,1}' \
		> "${region_specific_input_points_file}"

		# Check if points were actually extracted for this region
		if [ ! -f "${region_specific_input_points_file}" ] || ! [ -s "${region_specific_input_points_file}" ]; then
			echo "No points found for region ${region} for ${prefix}. Skipping this region."
			# Optionally remove the empty file if awk created one: rm -f "${region_specific_input_points_file}"
			continue # Skip to the next region in TARGET_REGIONS
		fi
		echo "Created ${region_specific_input_points_file} with point data."

		# --- 2. Adapt your first IF block for the current region ---
		# This creates the "modified" points file for 3dUndump input.
		# Output: space-separated i j k value (transformed coordinates)
		region_specific_modified_points_file="${folder}/new_points_modified_for_${region}.csv"
		#region_specific_modified_points_file="${folder}/points_modified_for_${region}.csv"
		rm $region_specific_modified_points_file

		master_nifti_for_undump="${folder}/${prefix}_rmprage.nii" # Master NIFTI, should exist
		z_dim_half=$(( $(3dinfo -nk $master_nifti_for_undump) / 2 ))

		echo $z_dim_half

		# Condition from your script: Create modified file if it doesn't exist AND the input points file (now region-specific) exists.
		# The existence and non-emptiness of region_specific_input_points_file is already confirmed above.
		if [ ! -f "${region_specific_modified_points_file}" ]; then
			# Awk transformation:
			# Input (region_specific_input_points_file) columns: $1=x, $2=y, $3=z, $4=t (comma-separated)
			# Transformed output for 3dUndump -ijk: i j k value --> x (216-y) (216-z) t (space-separated)
			#awk -F, '{printf "%s %s %s %s\n", (216-$2), (216-$1), $3, $4}' "${region_specific_input_points_file}" > "${region_specific_modified_points_file}"
			
			#echo awk -F, '{printf "%s %s %s %s\n", (216-$2), (216-$1), $z_dim_half, $4}' "${region_specific_input_points_file}" > "${region_specific_modified_points_file}"
			awk -v zval="$z_dim_half" -F, 'BEGIN{OFS=" "} {printf "%s %s %s %s\n", (216-$2), (216-$1), zval, $4}' "${region_specific_input_points_file}" > "${region_specific_modified_points_file}"

			echo "Created ${region_specific_modified_points_file}"
		else
			echo "${region_specific_modified_points_file} already exists. Skipping creation."
		fi

	done # End of loop through TARGET_REGIONS
done

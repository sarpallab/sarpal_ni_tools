template="/ihome/crc/install/fsl/6.0.4/centos7/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
#left_mask="left_hemisphere_mask.nii.gz"
#right_mask="right_hemisphere_mask.nii.gz"

#for folder in /ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/sub-2346/ses-20181221/; do
for folder in /ix1/ginger/dsarpal/lab/reorg/projects/mrsi/05_mrsi/SCZ/sub*/ses*/; do

	ref=$(basename $template)

        ln -sf ${template} ${folder}/${ref}
	
	sub=$(basename $(dirname $folder))
	ses=$(basename $folder)
	prefix=${sub}_${ses}

	scout_prefix=${prefix}"_scout_tonative_"
	native_prefix=${prefix}"_native_tomni_"
	
	native_brain="/ix1/ginger/dsarpal/lab/reorg/projects/7T/bids_correct/output/${sub}/${ses}/mprage/mprage_brain.nii.gz"

        ln -sf ${native_brain} ${folder}/${prefix}_mprage_brain.nii.gz

	if [ ! -f ${folder}/${prefix}_mprage.nii ] && [ -f ${folder}/mprage.nii ]; then
		${folder}/mprage.nii ${folder}/${prefix}_mprage.nii
	fi

	if [ ! -f ${folder}/${prefix}_rmprage.nii ] && [ -f ${folder}/rmprage.nii ]; then
		${folder}/rmprage.nii ${folder}/${prefix}_rmprage.nii
	fi

        if [ ! -f ${folder}/${scout_prefix}Warped.nii.gz ] &&  [ -f ${folder}/${prefix}_rmprage.nii ]; then
                echo antsRegistrationSyNQuick.sh -d 3 -m ${folder}/${prefix}_rmprage.nii -f ${folder}/${prefix}_mprage.nii -o ${folder}/${scout_prefix}
                antsRegistrationSyNQuick.sh -d 3 -m ${folder}/${prefix}_rmprage.nii -f ${folder}/${prefix}_mprage.nii -o ${folder}/${scout_prefix}
        fi

        if [ ! -f ${folder}/${native_prefix}Warped.nii.gz ] &&  [ -f ${folder}/${prefix}_rmprage.nii ]; then
                echo antsRegistrationSyNQuick.sh -d 3 -m ${folder}/${prefix}_mprage.nii -f ${folder}/${ref} -o ${folder}/${native_prefix}
                antsRegistrationSyNQuick.sh -d 3 -m ${native_brain} -f ${folder}/${ref} -o ${folder}/${native_prefix}
        fi

:<<'com'
	if [ ! -f ${folder}/points_modified.csv ] && [ -f ${folder}/points.csv ]; then
		awk -F, '{print (216-$3)" "(216-$2)" "$4" "$6}' ${folder}/points.csv > ${folder}/points_modified.csv
	fi

	if [ ! -f ${folder}/${prefix}_undumppoints.nii.gz ] && [ -f ${folder}/points.csv ]; then
		#3dUndump -prefix ${folder}/${prefix}_undumppoints.nii.gz -master ${folder}/${prefix}_rmprage.nii -srad 10 -cubes -overwrite -ijk <(sed 1d ${folder}/points_modified.csv|cut -d, -f1-3|paste - <(seq 1 2))
		3dUndump -prefix ${folder}/${prefix}_undumppoints.nii.gz -master ${folder}/${prefix}_rmprage.nii -srad 1 -cubes -overwrite -ijk <(tr -d '\r' < ${folder}/points_modified.csv)
		#tr -d '\r' < ${folder}/points_modified.csv
	fi
com

	TARGET_REGIONS=("lCaud" "rCaud" "lThal" "rThal")

	for region in "${TARGET_REGIONS[@]}"; do
    
		echo
		echo "Processing region: ${region} for ${prefix}"

		# --- 1. Prepare region-specific input point file ---
		# This file will contain original x,y,z,t from the master file for the current SUBJ, SES, and ROI.
		MASTER_POINTS_FILE_PATH="./points.csv"
		echo $MASTER_POINTS_FILE_PATH
		# Columns in this file will be: x,y,z,t (comma-separated)
		region_specific_input_points_file="${folder}/points_input_for_${region}.csv"
		echo $region_specific_input_points_file
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
		region_specific_modified_points_file="${folder}/next_points_modified_for_${region}.csv"
		#region_specific_modified_points_file="${folder}/points_modified_for_${region}.csv"
		rm $region_specific_modified_points_file

		master_nifti_for_undump="${folder}/${prefix}_rmprage.nii" # Master NIFTI, should exist
		z_dim_half=$(( $(3dinfo -nk $master_nifti_for_undump) / 2 ))

		# Condition from your script: Create modified file if it doesn't exist AND the input points file (now region-specific) exists.
		# The existence and non-emptiness of region_specific_input_points_file is already confirmed above.
		if [ ! -f "${region_specific_modified_points_file}" ]; then
			# Awk transformation:
			# Input (region_specific_input_points_file) columns: $1=x, $2=y, $3=z, $4=t (comma-separated)
			# Transformed output for 3dUndump -ijk: i j k value --> x (216-y) (216-z) t (space-separated)
			#awk -F, '{printf "%s %s %s %s\n", (216-$2), (216-$1), $3, $4}' "${region_specific_input_points_file}" > "${region_specific_modified_points_file}"
			awk -v zval="$z_dim_half" -F, 'BEGIN{OFS=" "} {printf "%s %s %s %s\n", (216-$2), (216-$1), zval, $4}' "${region_specific_input_points_file}" > "${region_specific_modified_points_file}"

			echo "Created ${region_specific_modified_points_file}"
		else
			echo "${region_specific_modified_points_file} already exists. Skipping creation."
		fi

		# --- 3. Adapt your second IF block for the current region ---
		# This runs 3dUndump using the "modified" points file.
		#region_specific_undumppoints_nii="${folder}/${prefix}_${region}_undumppoints.nii.gz"
		region_specific_undumppoints_nii="${folder}/next_${prefix}_${region}_undumppoints.nii.gz"
		master_nifti_for_undump="${folder}/${prefix}_rmprage.nii" # Master NIFTI, should exist

		# Condition from your script: Run 3dUndump if output NIFTI doesn't exist AND the original points file (now region-specific input) existed.
		# We also need to ensure the modified points file exists and is not empty.
		if [ -f "${region_specific_modified_points_file}" ] && [ -s "${region_specific_modified_points_file}" ] && [ ! -f "${region_specific_undumppoints_nii}" ]; then
			# Original commented line for 3dUndump using sed/cut/paste:
			#3dUndump -prefix ... -ijk <(sed 1d ${folder}/points_modified.csv|cut -d, -f1-3|paste - <(seq 1 2))
			# Original active line for 3dUndump:
			# 3dUndump -prefix ... -ijk <(tr -d '\r' < ${folder}/points_modified.csv)

			# Using the current region-specific modified points file
			3dUndump -prefix "${region_specific_undumppoints_nii}" \
				-master "${folder}/${prefix}_rmprage.nii" \
				-srad 2 -cubes -overwrite \
				-ijk <(tr -d '\r' < "${region_specific_modified_points_file}")
            
			echo "Created ${region_specific_undumppoints_nii}"
            
			# Original commented tr line (seems redundant if used in process substitution):
			# #tr -d '\r' < ${folder}/points_modified.csv 

		elif [ -f "${region_specific_undumppoints_nii}" ]; then
			echo "${region_specific_undumppoints_nii} already exists. Skipping 3dUndump for region ${region}."
		elif ! [ -s "${region_specific_modified_points_file}" ] || ! [ -f "${region_specific_modified_points_file}" ] ; then
			echo "Skipping 3dUndump for region ${region}: ${region_specific_modified_points_file} is empty or not found."
		fi
		echo "Finished processing for region: ${region} for ${prefix}"

	done # End of loop through TARGET_REGIONS


	TARGET_REGIONS=("lCaud" "rCaud" "lThal" "rThal")

	for region in "${TARGET_REGIONS[@]}"; do
		#region_specific_undumppoints="${folder}/${prefix}_${region}_undumppoints.nii.gz"
		region_specific_undumppoints="${folder}/next_${prefix}_${region}_undumppoints.nii.gz"
		#region_specific_undumppoints_MNI="${folder}/${prefix}_${region}_MNI.nii.gz"
		region_specific_undumppoints_MNI="${folder}/next_${prefix}_${region}_MNI.nii.gz"

		if [ ! -f ${region_specific_undumppoints_MNI} ] && [ -f ${region_specific_undumppoints} ]; then
			echo "Warping: ${region} for ${prefix}"
			#3dUndump -prefix ${folder}/${prefix}_undumppoints_mni.nii.gz -master ${folder}/${ref} -srad 10 -cubes -overwrite -xyz <(sed 1d ${folder}/points_mni.csv|cut -d, -f1-3|paste - <(seq 1 2))
                	antsApplyTransforms -d 3 \
				-i ${region_specific_undumppoints} \
				-o ${region_specific_undumppoints_MNI} \
				-r ${folder}/${ref} -n NearestNeighbor \
				-t ${folder}/${scout_prefix}1Warp.nii.gz \
				-t ${folder}/${scout_prefix}0GenericAffine.mat \
				-t ${folder}/${native_prefix}1Warp.nii.gz \
				-t ${folder}/${native_prefix}0GenericAffine.mat

			3drefit -space MNI ${region_specific_undumppoints_MNI}
		fi
	done
:<<'com'

	if [ -f ${folder}/${prefix}_undumppoints_mni.nii.gz ]; then
		echo ${prefix}
		echo "lCaud"
     		3dCM -mask $left_mask ${folder}/${prefix}_undumppoints_mni.nii.gz 
		echo "rCaud"
     		3dCM -mask $right_mask ${folder}/${prefix}_undumppoints_mni.nii.gz 

		echo "lThal" >> combined.txt
     		3dCM -mask $left_mask ${folder}/${prefix}_undumppoints_mni.nii.gz >> combined.txt
		echo "rThal" >> combined.txt
     		3dCM -mask $left_mask ${folder}/${prefix}_undumppoints_mni.nii.gz >> combined.txt

		echo "" >> combined.txt
com
done

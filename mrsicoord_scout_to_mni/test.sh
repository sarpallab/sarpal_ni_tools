template="/ihome/crc/install/fsl/6.0.4/centos7/fsl/data/standard/MNI152_T1_1mm.nii.gz"
left_mask="left_hemisphere_mask.nii.gz"
right_mask="right_hemisphere_mask.nii.gz"

#for folder in /ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/sub-2346/ses-20181221/; do
for folder in /ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/sub*/ses*/; do

	ref=$(basename $template)

        ln -sf ${template} ${folder}/${ref}

	scout_prefix="scout_to_native_"
	native_prefix="native_to_mni_"
	
	sub=$(basename $(dirname $folder))
	ses=$(basename $folder)
	prefix=${sub}_${ses}

	if [ ! -f ${folder}/${prefix}_mprage.nii ] && [ -f ${folder}/mprage.nii ]; then
		mv ${folder}/mprage.nii ${folder}/${prefix}_mprage.nii
	fi

	if [ ! -f ${folder}/${prefix}_rmprage.nii ] && [ -f ${folder}/rmprage.nii ]; then
		mv ${folder}/rmprage.nii ${folder}/${prefix}_rmprage.nii
	fi

        if [ ! -f ${folder}/${scout_prefix}1Warp.nii.gz ] &&  [ -f ${folder}/${prefix}_rmprage.nii ]; then
                antsRegistrationSyNQuick.sh -d 3 -m ${folder}/${prefix}_rmprage.nii -f ${folder}/${prefix}_mprage.nii -o ${folder}/${scout_prefix}
        fi

        if [ ! -f ${folder}/${native_prefix}1Warp.nii.gz ] &&  [ -f ${folder}/${prefix}_rmprage.nii ]; then
                antsRegistrationSyNQuick.sh -d 3 -m ${folder}/${prefix}_mprage.nii -f ${folder}/${ref} -o ${folder}/${native_prefix}
        fi

	if [ ! -f ${folder}/points_modified.csv ] && [ -f ${folder}/points.csv ]; then
		awk -F, '{print (216-$3)" "(216-$2)" "$4" "$6}' ${folder}/points.csv > ${folder}/points_modified.csv
	fi

	if [ ! -f ${folder}/${prefix}_undumppoints.nii.gz ] && [ -f ${folder}/points.csv ]; then
		#3dUndump -prefix ${folder}/${prefix}_undumppoints.nii.gz -master ${folder}/${prefix}_rmprage.nii -srad 10 -cubes -overwrite -ijk <(sed 1d ${folder}/points_modified.csv|cut -d, -f1-3|paste - <(seq 1 2))
		3dUndump -prefix ${folder}/${prefix}_undumppoints.nii.gz -master ${folder}/${prefix}_rmprage.nii -srad 1 -cubes -overwrite -ijk <(tr -d '\r' < ${folder}/points_modified.csv)
		#tr -d '\r' < ${folder}/points_modified.csv
	fi

:<<'EOF'
	if [ ${folder}/${prefix}_undumppoints.nii.gz ]; then
		for opt in ""; do
     			3dCM -roi_vals 1 2 $opt ${folder}/${prefix}_undumppoints.nii.gz 2>/dev/null|sed 1d|paste - - 
     			#3dCM -roi_vals 1 2 $opt ${folder}/${prefix}_undumppoints.nii.gz 2>/dev/null|sed 1d|paste | awk '{print $1 "," $2 "," $3 ",0," NR}' > ${folder}/temp_points.csv 
			#3dCM -roi_vals 1 2 $opt ${input_csv} 2>/dev/null | sed '/^#/d' | sed 1d | awk '{print $1 "," $2 "," $3 ",0," NR % 2 == 1 ? 1 : 2}' > ${folder}/temp_points.csv 
			#3dCM -roi_vals 1 2 $opt ${input_csv} 2>/dev/null | sed '/^#/d' | sed 1d | awk '{print $1 "," $2 "," $3 ",0," NR}' > ${folder}/temp_points.csv
			#3dCM -roi_vals 1 2 $opt ${input_csv} 2>/dev/null | sed 1d | awk '!/^#/{print $1 "," $2 "," $3 ",0," NR}' > ${folder}/temp_points.csv
   		done
	fi

	#output_csv=xyz_points.csv

	#echo "x,y,z,t,label" > ${folder}/${output_csv}	
	#cat ${folder}/temp_points.csv >> ${folder}/${output_csv}

        if [ ! -f ${folder}/points_mni.csv ]; then

                antsApplyTransformsToPoints -d 3 \
			-i ${folder}/xyz_points.csv \
			-o ${folder}/points_mni.csv \
			-t ${folder}/${scout_prefix}1Warp.nii.gz \
			-t ${folder}/${scout_prefix}0GenericAffine.mat \
			-t ${folder}/${native_prefix}1Warp.nii.gz \
			-t ${folder}/${native_prefix}0GenericAffine.mat
        fi
EOF

	if [ ! -f ${folder}/${prefix}_undumppoints_mni.nii.gz ] && [ -f ${folder}/${prefix}_undumppoints.nii.gz ]; then
		#3dUndump -prefix ${folder}/${prefix}_undumppoints_mni.nii.gz -master ${folder}/${ref} -srad 10 -cubes -overwrite -xyz <(sed 1d ${folder}/points_mni.csv|cut -d, -f1-3|paste - <(seq 1 2))
                antsApplyTransforms -d 3 \
			-i ${folder}/${prefix}_undumppoints.nii.gz \
			-o ${folder}/${prefix}_undumppoints_mni.nii.gz \
			-r ${folder}/${ref} -n NearestNeighbor \
			-t ${folder}/${scout_prefix}1Warp.nii.gz \
			-t ${folder}/${scout_prefix}0GenericAffine.mat \
			-t ${folder}/${native_prefix}1Warp.nii.gz \
			-t ${folder}/${native_prefix}0GenericAffine.mat

		3drefit -space MNI ${folder}/${prefix}_undumppoints_mni.nii.gz
	fi

	if [ -f ${folder}/${prefix}_undumppoints_mni.nii.gz ]; then
		echo ${prefix}
		echo "lSTS"
     		3dCM -mask $left_mask ${folder}/${prefix}_undumppoints_mni.nii.gz 
		echo "rSTS"
     		3dCM -mask $right_mask ${folder}/${prefix}_undumppoints_mni.nii.gz 

		echo ${prefix} >> combined.txt
		echo "lSTS" >> combined.txt
     		3dCM -mask $left_mask ${folder}/${prefix}_undumppoints_mni.nii.gz >> combined.txt

		echo "" >> combined.txt
	fi

done

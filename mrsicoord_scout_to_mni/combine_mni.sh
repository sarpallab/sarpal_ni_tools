template="/ihome/crc/install/fsl/6.0.4/centos7/fsl/data/standard/MNI152_T1_1mm.nii.gz"

#for folder in /ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/sub-2346/ses-20181221/; do
#for folder in /ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/sub-2351/ses-20190124/; do

session_files=()

for folder in /ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/sub*/ses*/; do

	ref=$(basename $template)
		
	first_ses=$(find "$folder" -type d -name "ses-*" | sort  | head -n 1)

	if [ -n "$first_ses" ]; then

		# Find the NIFTI file within the first session folder
		nifti_file=$(find $first_ses -type f -name "*undumppoints_mni.nii.gz" | head -n 1)
    
		if [ -f $nifti_file ]; then
			session_files+=("$nifti_file")
		fi
	fi
done

echo ${session_files[@]}

outfile="/ix1/ginger/dsarpal/lab/data/7T/mrsi_reorg/result.nii.gz"
# Initialize the output file by copying the first session file of the first subject
cp ${session_files[0]} $outfile

# Add each subsequent NIFTI file to the output file
for nifti_file in ${session_files[@]:1}; do
	fslmaths $outfile -add $nifti_file $outfile
done

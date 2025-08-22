import pandas as pd
import os
import subprocess

def clusterize(funcfile, invertedfile, prefix, native_target):
	
	invert_command=f"fslmaths {funcfile} -mul -1 {invertedfile}.nii.gz"
	
	if not os.path.exists(f"{invertedfile}.nii.gz"):
		subprocess.run(invert_command, shell=True)

	
	mask_command=f"fslmaths {invertedfile}.nii.gz -mas {native_target} {invertedfile}_masked.nii.gz"

	if not os.path.exists(f"{invertedfile}_masked.nii.gz"):
		subprocess.run(mask_command, shell=True)


	cluster_command=f"cluster -i {invertedfile}_masked.nii.gz --thresh=0.30 --oindex={prefix}index.nii.gz --omean={prefix}mean.nii.gz --osize={prefix}size.nii.gz > {prefix}info.txt"
	print(cluster_command)

	if not os.path.exists(f"{prefix}index.nii.gz"):
		subprocess.run(cluster_command, shell=True)

	
	

subj="sub-1019"
analysis_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/bids_mni/derivatives/analysis_mni/{subj}"

funcfile=f"{analysis_dir}/{subj}_basalis_corrmap_z.nii.gz"
invertedfile=f"{analysis_dir}/{subj}_basalis_corrmap_z_inverse"
native_dlpfc=f"{analysis_dir}/{subj}_funcmni_dlpfc.nii.gz"
prefix=f"{analysis_dir}/{subj}_basalis_cluster_"

#funcfile=f"{analysis_dir}/{subj}_subgenual_corrmap_z.nii.gz"
#invertedfile=f"{analysis_dir}/{subj}_subgenual_corrmap_z_inverse"
#native_dlpfc=f"{analysis_dir}/{subj}_funcmni_dlpfc.nii.gz"
#prefix=f"{analysis_dir}/{subj}_subgenual_cluster_"

#clusterize(funcfile, invertedfile, prefix, native_dlpfc)
clusterize(funcfile, invertedfile, prefix, native_dlpfc)

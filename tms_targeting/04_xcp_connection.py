import pandas as pd
import os
import subprocess

def move_files(subj,xcp_dir, analysis_dir):

	dirdict = {'anat':'space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz',
			'func':'task-rest_space-MNI152NLin2009cAsym_desc-denoisedSmoothed_bold.nii.gz' }
	
	os.makedirs(analysis_dir, exist_ok=True)
	
	for dirname, filename in dirdict.items():
		old_name = f"{xcp_dir}/{dirname}/{subj}_{filename}"
		new_name = f"{analysis_dir}/{subj}_{filename}"
		if not os.path.exists(new_name):
			subprocess.run(f"ln -sf {old_name} {new_name}", shell=True)
	

def get_roi(roi, func, destination):
	
	#basalis_func=f"{analysis_dir}/{subj}_funcmni_fnirt_basalis.nii.gz"
	#resample_command = f"flirt -in {roi} -ref {func} -out {basalis_func}"

	resample_command = f"3dresample -input {roi} -master {func} -prefix {destination}"
	if not os.path.exists(destination):
		subprocess.run(resample_command,shell=True)
	



def get_connection_map(prefix, funcfile, anatfile, region, region_name):
	
	prefix += f"_{region_name}"
	oneD = f"{prefix}.1D"

	maskave_command=f"3dmaskave -q -mask {region} {funcfile} > {oneD}"
	if not os.path.exists(oneD):
		subprocess.run(maskave_command, shell=True)
		#subprocess.run(fslmaths_command, shell=True)
		#subprocess.run(fslmeants_command, shell=True)

	
	corrmap_r = f"{prefix}_corrmap_r.nii.gz"
	Tcorr_command=f"fsl_glm -d {oneD} -i {funcfile} -o {corrmap_r}"
	if not os.path.exists(corrmap_r):
		subprocess.run(Tcorr_command, shell=True)
	


	corrmap_z = f"{prefix}_corrmap_z.nii.gz"
	corrz_command=f"3dcalc -a {corrmap_r} -expr 'log((1+a)/(1-a))/2' -prefix {corrmap_z}"
	if not os.path.exists(corrmap_z):
		subprocess.run(corrz_command, shell=True)


def smooth_connection(prefix, funcfile, region_name):

	prefix += f"_{region_name}"
	corrmap_z = f"{prefix}_corrmap_z.nii.gz"
	outfile = f"{prefix}_corrmap_z_smoothed.nii.gz"
	
	prefix += f"_{region_name}"
	smooth_command = f"3dBlurToFWHM -FWHM 8 -prefix {outfile} -input {corrmap_z}"

	if not os.path.exists(outfile):
		subprocess.run(smooth_command, shell=True)

subj="sub-1019"

xcp_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/bids_mni/xcp_d/{subj}"
analysis_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/bids_mni/derivatives/analysis_mni/{subj}"
roi_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/scripts/rois/"

basalis=f"{roi_dir}/basalis_1-4.nii.gz"
subgenual=f"{roi_dir}/subgenual.nii.gz"
dlpfc=f"{roi_dir}/dlpfc.nii"

funcfile=f"{analysis_dir}/{subj}_task-rest_space-MNI152NLin2009cAsym_desc-denoisedSmoothed_bold.nii.gz"
anatfile=f"{analysis_dir}/{subj}_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz"
prefix=f"{analysis_dir}/{subj}"

dlpfc_func=f"{analysis_dir}/{subj}_funcmni_dlpfc.nii.gz"
basalis_func=f"{analysis_dir}/{subj}_funcmni_basalis.nii.gz"
subgenual_func=f"{analysis_dir}/{subj}_funcmni_subgenual.nii.gz"

#Check if exists
move_files(subj, xcp_dir, analysis_dir)


get_roi(dlpfc,funcfile,dlpfc_func)
get_roi(basalis,funcfile,basalis_func)
get_roi(subgenual,funcfile,subgenual_func)

#get_connection_map(prefix, funcfile, basalis_func, "basalis")
get_connection_map(prefix, funcfile, anatfile, basalis_func, "basalis")
get_connection_map(prefix, funcfile, anatfile, subgenual_func, "subgenual")

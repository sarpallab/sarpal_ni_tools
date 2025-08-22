import pandas as pd
import os
import subprocess

'''
def prepare_warp_funcavg(analysis_dir, native_anat, func_file, avg_file, out_prefix):
	
	if not os.path.exists(avg_file):
		average_command = f"fslmaths {func_file} -Tmean {avg_file}"
		subprocess.run(average_command, shell=True)

	if not os.path.exists(f"{out_prefix}1warp.nii.gz"):
		register_command=f"./antsRegistrationSyNQuick.sh -d 3 -f {native_anat} -m {avg_file} -o {out_prefix} -t s"
		subprocess.run(register_command, shell=True)
'''

def prepare_warp(analysis_dir, native_anat, mnianat_file, out_prefix):
	

	if not os.path.exists(f"{out_prefix}1warp.nii.gz"):
		register_command=f"antsRegistrationSyNQuick.sh -d 3 -f {native_anat} -m {mnianat_file} -o {out_prefix} -t s"
		subprocess.run(register_command, shell=True)


def apply_warp(analysis_dir, subj, ref_file, to_warp, out_prefix, outfile):

	if not os.path.exists(f"{outfile}.nii.gz"):
		warp_command=f"antsApplyTransforms --verbose -d 3 -n NearestNeighbor -i {to_warp} -r {ref_file} -t {out_prefix}1Warp.nii.gz -t {out_prefix}0GenericAffine.mat -o {outfile}.nii.gz"
		subprocess.run(warp_command, shell=True)

	binarize_command=f"fslmaths {outfile}.nii.gz -bin {outfile}.nii.gz"
	subprocess.run(binarize_command, shell=True)

	#resample_command=f"3dresample -prefix {outfile}_res.nii.gz -input {outfile}.nii.gz -master {avg_file}"
	#subprocess.run(resample_command, shell=True)





subj="sub-1020"

fmri_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/bids_mni/derivatives/{subj}"
xcp_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/bids_mni/derivatives/{subj}"
analysis_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/bids_mni/derivatives/analysis_mni/{subj}"
out_dir=f"/ix1/ginger/dsarpal/lab/reorg/projects/enact/bids_mni/derivatives/out/{subj}"
roi_dir="/ix1/ginger/dsarpal/lab/reorg/projects/enact/scripts/rois"

#nativeanatfile=f"{fmri_dir}/anat/{subj}_desc-brain_mask.nii.gz"
#newnativeanatfile=f"{out_dir}/{subj}_desc-brain_mask.nii.gz" 

nativeanatfile=f"{fmri_dir}/anat/{subj}_desc-preproc_T1w.nii.gz"
newnativeanatfile=f"{out_dir}/{subj}_desc-preproc_T1w.nii.gz" 
#funcfile=f"{analysis_dir}/{subj}_task-rest_space-MNI152NLin2009cAsym_desc-denoisedSmoothed_bold.nii.gz"
#newfuncfile=f"{out_dir}/{subj}_task-rest_space-MNI152NLin2009cAsym_desc-denoisedSmoothed_bold.nii.gz"
mnianatfile=f"{fmri_dir}/anat/{subj}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
newmnianatfile=f"{out_dir}/{subj}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"

anatbrainfile=f"{fmri_dir}/anat/{subj}_desc-brain_mask.nii.gz"
newanatbrainfile=f"{out_dir}/{subj}_desc-brain_mask.nii.gz"

#funcbrainfile=f"{analysis_dir}/{subj}_task-rest_space-MNI152NLin2009cAsym_desc-denoisedSmoothed_bold.nii.gz"
#newfuncbrainfile=f"{out_dir}/{subj}_task-rest_space-MNI152NLin2009cAsym_desc-denoisedSmoothed_bold.nii.gz"


if not os.path.islink(newnativeanatfile): os.symlink(nativeanatfile, newnativeanatfile)
if not os.path.islink(newmnianatfile): os.symlink(mnianatfile, newmnianatfile)
if not os.path.islink(newanatbrainfile): os.symlink(anatbrainfile, newanatbrainfile)
#if not os.path.islink(newfuncfile): os.symlink(funcfile, newfuncfile)


avg_file = f"{out_dir}/{subj}_func_avg.nii.gz"
#outfile=f"{analysis_dir}/{subj}/{subj}_task-rest_space-T1w_desc-preproc_bold.nii.gz"

mni_template=f"{roi_dir}/mni_icbm152_t1_tal_nlin_asym_09c_brain_2.3mm.nii"
#out_prefix=f"{out_dir}/{subj}_mni_to_T1_funcavg_"
out_prefix=f"{out_dir}/{subj}_mni_to_anat_"


cluster_to_warp=f"{analysis_dir}/cluster2.nii.gz"
clustermap_to_warp=f"{analysis_dir}/{subj}_basalis_corrmap_z.nii.gz"
cluster_native=f"{out_dir}/{subj}_cluster2"
clustermap_native=f"{out_dir}/{subj}_clustermap"

'''
dlpfc_to_warp=f"{analysis_dir}/{subj}_funcmni_dlpfc.nii"
subgenual_to_warp=f"{roi_dir}/{subj}_funcmni_subgenual.nii.gz"
basalis_to_warp=f"{roi_dir}/{subj}_funcmni_basalis.nii.gz"

dlpfc_native_funcavg=f"{out_dir}/{subj}_native_funcavg_dlpfc"
basalis_native_funcavg=f"{out_dir}/{subj}_native_funcavg_basalis"
subgenual_native_funcavg=f"{out_dir}/{subj}_native_funcavg_subgenual"
'''
if not os.path.exists(f"{out_prefix}1Warp.nii.gz"):
	os.makedirs(out_dir, exist_ok = True)
	#prepare_warp_funcavg(analysis_dir, anatfile, funcfile, avg_file, out_prefix)
	prepare_warp(analysis_dir, newanatbrainfile, newmnianatfile, out_prefix)

if not os.path.exists(f"{cluster_native}.nii.gz"):
	apply_warp(out_dir, subj, nativeanatfile, cluster_to_warp, out_prefix, cluster_native)
#if not os.path.exists(f"{clustermap_native}.nii.gz"):
	#apply_warp(out_dir, subj, anatfile, clustermap_to_warp, out_prefix, clustermap_native)

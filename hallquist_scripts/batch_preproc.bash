#!/usr/bin/env bash

#SBATCH -A dsarpal
#SBATCH --ntasks-per-node 1
#SBATCH --nodes 1
#SBATCH --time 3:00:00

# run for 3 hours. should finish under 2:30

# USAGE
#  command line
#    ./batch_preproc.bash /ix1/ginger/dsarpal/lab/reorg/projects/7T/bids_test/sub-2808/ses-20230818/func/sub-2808_ses-20230818_task-rest_bold.nii.gz 
# batch
#    sbatch -J 2808_20230818 --export ALL,INFILE=/ix1/ginger/dsarpal/lab/reorg/projects/7T/bids_test/sub-2808/ses-20230818/func/sub-2808_ses-20230818_task-rest_bold.nii.gz 

set -euo pipefail
# needed modules. NB. don't load  r/4.1.0. afni loads 4.0.0
<<'###'

module load afni fsl ants convert3d matlab
# make sure each node ica aroma
pip show icaaroma || 
  python3 -m pip install --user future git+https://github.com/WillForan/ICA-AROMA.git@maartenmennes-setup.py

# setup R: need writeable R_LIBS_USER dir and to make a CRAN repo (CMU) known
# all for orthopolynom, only needed for resting state bandpass+regression of csf, wm, motion, etc
export R_LIBS_USER=/ix1/ginger/dsarpal/ni_tools/R-4.0.0
mkdir -p "$R_LIBS_USER"
grep -q R_LIBS_USER ~/.Renviron || echo "R_LIBS_USER=$R_LIBS_USER" >> $_
grep -q CRAN ~/.Rprofile || echo "options(repos=structure(c(CRAN='http://lib.stat.cmu.edu/R/CRAN')))" >> $_
 #grep -q $R_LIBS_USER ~/.Rprofile || echo ".libPaths('$R_LIBS_USER')" >> $_;
Rscript -e 'if(!require("orthopolynom")) install.packages("orthopolynom")'
export PATH="/ix1/ginger/dsarpal/lab/ni_tools/ROBEX:/ix1/ginger/dsarpal/lab/reorg/ni_tools/fmri_processing_scripts:$PATH"
export MRI_STDDIR=/ix1/ginger/dsarpal/lab/ni_tools/standard/

# afnis 3dREMLfit fails to write data on CRC w/more than 1 thread?
# only have 1 core anyway
export OMP_NUM_THREADS=1
###

source '/ix1/ginger/dsarpal/lab/reorg/projects/toronto_processing/scripts/preproc_setup.bash'

# preprocessFunctional -check_dependencies # okay missing nipy (no 4dslice here) and R. but need aroma

 bidsdir_input='/ix1/ginger/dsarpal/lab/reorg/projects/toronto_processing/preprocess'
        outdir='/ix1/ginger/dsarpal/lab/reorg/projects/toronto_processing/preprocess'

# with global envir var input via slurm (or for testing, first argument in)
[ $# -eq 1 ] && INFILE="$1" || INFILE="${INFILE:-}" # set default to empty string to avoid 'set -u' error
[ -z "$INFILE" -o ! -s "$INFILE" ] && echo "ERROR no inputfile INFILE ('$INFILE')" >&2 && exit 1
INFILE="$(readlink -f "$INFILE")"
JSONFILE="$(echo $INFILE | sed 's/.nii/.json/')"

# find mpage dir and check input files exist
! [[ $INFILE =~ sub-([A-Z]{2})([0-9]{3})* ]]  && echo "ERROR INFILE ('$INFILE') doesn't have sub-*/*" >&2 && exit 1
idpath=$BASH_REMATCH
idses=${idpath}

# orig mprage
anat_dir="$bidsdir_input/$idpath/anat"
#echo "# trying bids run preproc in $mprage_dir"
warp_coef="$anat_dir/${idses}_T1w_warpcoef.nii.gz"
mprage_bet="$anat_dir/${idses}_T1w_bet.nii.gz"

echo $warp_coef

# 2023-11 WS version
if [ ! -r "$warp_coef" ]; then
  mprage_dir="$bidsdir_input/$idpath/anat"
  echo "# WARNING: cannot find '$warp_coef'; trying bids run preproc in $anat_dir"
  warp_coef="$anat_dir/${idses}_T1w_warpcoef.nii.gz"
  mprage_bet="$anat_dir/${idses}_T1w_bet.nii.gz"
fi

[ ! -r "$warp_coef" ] && echo "ERROR: cannot find warpcoef '$warp_coef'" >&2 && exit 2
[ ! -r "$mprage_bet" ] && echo "ERROR: cannot find mprage_bet '$mprage_bet'" >&2 && exit 2

# actually run
preproc_dir=$outdir/$idpath/rest/
mkdir -p $preproc_dir 
cd $preproc_dir

! test -r ${idses}_rest.nii.gz && ln -s $INFILE ${idses}_rest.nii.gz
#! test -r ${idses}_rest.json && ln -s $JSONFILE ${idses}_rest.json

export MATLAB_RAM_limit=2

echo y | preprocessFunctional \
  -4d ${idses}_rest.nii.gz \
  -mprage_bet "$mprage_bet" \
  -warpcoef "$warp_coef" \
  -template_brain MNI_2mm \
  -slice_acquisition seqdesc -ica_aroma \
  -st_first \
  -tr 2 \
  -wavelet_despike \
  -wavelet_threshold 10   -threshold 98_2 \
  -bandpass_filter 0.009   .08 \
  -rescaling_method       10000_globalmedian \
  -func_struc_dof bbr \
  -warp_interpolation     spline \
  -constrain_to_template  y \
  -motion_censor  fd=0.3,dvars=20 \
  -smoothing_kernel       5 \
  -nuisance_file  nuisance_regressors.txt \
  -nuisance_compute       csf,dcsf,wm,dwm \
  -nuisance_regression    6motion,d6motion,csf,dcsf,wm,dwm

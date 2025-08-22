#!/usr/bin/env bash

#SBATCH -A dsarpal
#SBATCH --ntasks-per-node 1
#SBATCH --nodes 1
#SBATCH --time 2:30:00

# run for 2 hours and 30 seconds

# USAGE
#  command line
#    ./batch_preproc_t1.bash /ix1/ginger/dsarpal/lab/reorg/projects/7T/bids_test/sub-2808/ses-20230818/anat/sub-2808_ses-20230818_T1w.nii.gz
# batch
#    sbatch -J 2808_20230818 --export ALL,INFILE=/ix1/ginger/dsarpal/lab/reorg/projects/7T/bids_test/sub-2808/ses-20230818/anat/sub-2808_ses-20230818_T1w.nii.gz ./batch_preproc_t1.bash 

set -euo pipefail

<<'###'
module load afni fsl ants convert3d 

export PATH="/ix1/ginger/dsarpal/lab/ni_tools/ROBEX:/ix1/ginger/dsarpal/lab/reorg/programs/fmri_processing_scripts:$PATH"
export MRI_STDDIR=/ix1/ginger/dsarpal/lab/ni_tools/standard/
preprocessMprage -check_dependencies
###


 source '/ix1/ginger/dsarpal/lab/reorg/projects/toronto_processing/scripts/preproc_setup.bash'

 bidsdir_input='/ix1/ginger/dsarpal/lab/reorg/projects/toronto_processing/toronto_data'
        outdir='/ix1/ginger/dsarpal/lab/reorg/projects/toronto_processing/preprocess'

# with global envir var input via slurm (or for testing, first argument in)
[ $# -eq 1 ] && INFILE="$1" || INFILE="${INFILE:-}" # set default to empty string to avoid 'set -u' error
[ -z "$INFILE" -o ! -s "$INFILE" ] && echo "ERROR no inputfile INFILE ('$INFILE')" >&2 && exit 1
INFILE="$(readlink -f "$INFILE")"
echo $INFILE

# find mpage dir and check input files exist
! [[ $INFILE =~ sub-([A-Z]{2})([0-9]{3})* ]]  && echo "ERROR INFILE ('$INFILE') doesn't have sub-*/ses*" >&2 && exit 1
idpath=$BASH_REMATCH
idses=${idpath}

echo $idses

# actually run
preproc_dir=$outdir/$idpath/anat/
mkdir -p $preproc_dir 
cd $preproc_dir
! test -r ${idses}_T1w.nii.gz && cp $INFILE ./${idses}_T1w.nii.gz

preprocessMprage -n ${idses}_T1w.nii.gz

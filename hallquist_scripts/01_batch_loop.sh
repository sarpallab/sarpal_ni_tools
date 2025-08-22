#!/usr/bin/env bash
#
# run hallquist pipeline w/sbatch (mprage and rest reprocessing)
# only for not complete and not running
# mprage launched and rest made a dependency on that job
# see batch_preproc.bash and batch_preproc_t1.bash
#
# USAGE:
#  # only 2023 visits. only print what would submit
#  DRYRUN=1 ./01_batch_loop.sh /ix1/ginger/dsarpal/lab/reorg/projects/7T/bids_test/sub*/ses-2023*/func/*_task-rest_bold.nii.gz
#
#  # all for real
#  ./01_batch_loop.sh 

# 20231116WS - init
# 20231116WF - mprage dependencies. dryrun. account by user

set -euo pipefail

# run on what we're given or on everything
[ $# -gt 0 ] &&
  RESTFILES=("$@") ||
  RESTFILES=(/ix1/ginger/dsarpal/lab/reorg/projects/toronto_processing/toronto_data/sub*/func/*_task-rest_bold.nii)

# prefix command like 'DRYRUN=1 ./01_batch_loop.sh ...' to echo instead of actually running
[ ! -v DRYRUN ] && DRYRUN= || DRYRUN=echo
# TODO: remove echo when ready
#DRYRUN=echo

# foran cant submit to dsarpal
[ $USER == warren ] && slurmaccount='-A npac' || slurmaccount=""

mkdir -p log
for f in ${RESTFILES[@]}; do
  jobname=$(basename "$f" _task-rest_bold.nii.gz)

  # dont need to submit if already finished
  finalout=$(dirname $f|sed s/toronto_data/preprocess/)/.preprocessfunctional_complete
  test -r $finalout && echo "# already finished $jobname (have $finalout)" && continue

  # make sure we're not already running the job 
  squeue -o %j | grep -e "$jobname$" && echo "# already running $jobname (working on $finalout)" && continue

  # make batch_preproc dependent on batch_preproc_t1
  warpcoef=$(dirname $finalout|sed 's/toronto_data/preprocess/')/mprage_warpcoef.nii.gz
  deps=
  if [ ! -r $warpcoef ]; then

    anat_name=$(echo ${jobname} | sed 's/task-rest_bold/T1w/')
    mprage=$(dirname $f |sed s/func/anat/)/${anat_name}
    ! test -r $mprage && echo "ERROR: no t1 '$mprage'" && continue
    jobnamemprage=$jobname-mprage
    if ! squeue -o %j | grep "$jobnamemprage"; then
      echo "# queue for mprage: $warpcoef"
      $DRYRUN sbatch $slurmaccount -o log/%x-%u-%j.log \
                  -J "$jobnamemprage" \
                  --export=All,INFILE=$mprage \
                  ./batch_preproc_t1.bash
    fi
    mjid=$(squeue  -o "%j %i"|grep -Po "(?<=$jobnamemprage ).*"||:)
    # 'afterany' should maybe be 'afterok' -- but not sure if preprocessMprage exist status is always 0
    [ -n $mjid ] && deps="--dependency=afterany:$mjid"
  fi

  echo "# queue for rest: $finalout"
  $DRYRUN sbatch -o log/%x-%u-%j.log \
              -J "$jobname" \
              --nodes=1 --ntasks-per-node=8 --mem=16G \
              --export=ALL,INFILE=$f \
              $deps \
              ./batch_preproc.bash

  # only care to try/see one visit? break afer first iteration
done

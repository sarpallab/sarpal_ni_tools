#if [ -z "${CONDA_PYTHON_EXE:-}" ]; then
module load python/ondemand-jupyter-python3.10 ants afni singularity fsl
source activate /ix1/ginger/dsarpal/wcs17/envs/enact
#fi

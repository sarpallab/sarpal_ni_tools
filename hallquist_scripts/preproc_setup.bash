# 20231127WF - CRC setup for preprocessFunctional and preprocessMprage
#  use in scripts like
#    source /ix1/ginger/dsarpal/lab/reorg/projects/7T/scripts/preproc_setup.bash
#
#  check with
#   preprocessFunctional -check_dependencies
#   preprocessMprage     -check_dependencies
#

#module load vim/8.1
module load afni fsl ants convert3d matlab squashfs-tools libjpeg imagemagick
module load python/ondemand-jupyter-python3.10

python3 -m pip install nipy --upgrade --user

pip show icaaroma || 
  python3 -m pip install --user \
      future \
      git+https://github.com/WillForan/ICA-AROMA.git@maartenmennes-setup.py

# setup R: need writeable R_LIBS_USER dir and to make a CRAN repo (CMU) known
# all for orthopolynom, only needed for resting state bandpass+regression of csf, wm, motion, etc
export R_LIBS_USER=/ix1/ginger/dsarpal/ni_tools/R-4.0.0
mkdir -p "$R_LIBS_USER"
# THIS WILL BE WEIRD IF WE UPDATE R. only writes once. will still be 4.0.0
grep -q R_LIBS_USER ~/.Renviron || echo "R_LIBS_USER=$R_LIBS_USER" >> $_
grep -q CRAN ~/.Rprofile || echo "options(repos=structure(c(CRAN='http://lib.stat.cmu.edu/R/CRAN')))" >> $_

#grep -q $R_LIBS_USER ~/.Rprofile || echo ".libPaths('$R_LIBS_USER')" >> $_;
Rscript -e 'if(!require("orthopolynom")) install.packages("orthopolynom")'

export PATH="/ix1/ginger/dsarpal/lab/ni_tools/ROBEX:/ix1/ginger/dsarpal/lab/reorg/ni_tools/fmri_processing_scripts:$PATH:$HOME/.local/bin:/ix1/ginger/dsarpal/lab/reorg/programs/lncdtools"
export MRI_STDDIR=/ix1/ginger/dsarpal/lab/ni_tools/standard/

# only use a single core. e.g. afni's 3dREMLfit
export OMP_NUM_THREADS=1

[[ $(R --version|sed 1q|grep -Po '\d+\.\d+\.\d+') == 4.0.0 ]] || echo "WARNING: wrong R version!?"

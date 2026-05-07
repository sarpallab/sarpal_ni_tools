# 2026-04-28 AndyP
# NDA Script

library(tidyverse)

cols2add <- c("transmit_coil","image_history","qc_outcome","qc_description","qc_fail_quest_reason","decay_correction","frame_end_times","frame_end_unit","frame_start_times","frame_start_unit","pet_isotope","pet_tracer","time_diff_inject_to_image","time_diff_units","pulse_seq","slice_acquisition","software_preproc","study","week","experiment_description","year_mta","timepoint_label","aqi","fd_mean","dvars_std","tsnr","fetal_age","fetal_age_type","accession_number","ageyears","iti_onset","stim1","stim2","stim1_side","stim1_magnitude","stim2_magnitude","choice_side","computer_choice","stim1_outcome","session_fmri","choice_fmri","outcome_fmri","rt_fmri","task__version","block_sv","trial_num","options_onset","cue_onset","interval_onset","monitor_onset","session_det","gbc","vtca","vtcan","eventname","taskname","abbrev_taskname","vendor","gbc_r","gbc_region","pib_global","at_scanner","groups","mechanism_measure","session_id","image_extent5","extent5_type","image_unit5","image_resolution5","excitation_wavelength")
cols2add2 <- c("procdate","visnum","block_number","level","stain_details","comments_misc","image_thumbnail_file","transmit_coil","image_history","qc_outcome","qc_description","qc_fail_quest_reason","decay_correction","frame_end_times","frame_end_unit","frame_start_times","frame_start_unit","pet_isotope","pet_tracer","time_diff_inject_to_image","time_diff_units","pulse_seq","slice_acquisition","software_preproc","study","week","experiment_description","year_mta","timepoint_label","aqi","fd_mean","dvars_std","tsnr","fetal_age","fetal_age_type","accession_number","ageyears","iti_onset","stim1","stim2","stim1_side","stim1_magnitude","stim2_magnitude","choice_side","computer_choice","stim1_outcome","session_fmri","choice_fmri","outcome_fmri","rt_fmri","task__version","block_sv","trial_num","options_onset","cue_onset","interval_onset","monitor_onset","session_det","gbc","vtca","vtcan","eventname","taskname","abbrev_taskname","vendor","gbc_r","gbc_region","pib_global","at_scanner","groups","mechanism_measure","session_id")

df <- read.csv('/Users/andypapale/Downloads/image03.csv',skip=1)
df[cols2add] <- NA
df[cols2add2] <- NA
setwd('~/Desktop')
con <- file('image03.csv', open="wb")
writeLines("",con)
#write_csv(df,con,na = "",append=TRUE, col_names = TRUE)
write_csv(df, con,append=TRUE,col_names=TRUE)
close(con)

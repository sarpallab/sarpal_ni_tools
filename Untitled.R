# 2026-07-06 AndyP
# wrap_res_with_age.R
# wrapper function to get GAMs for MRSI



df <- read_csv('/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper/13MP20200207_LCMv2fixidx_Raw.csv')
df <- df %>% separate_wider_delim(cols = ld8,delim="_",names=c("id","dateNumeric"),cols_remove=FALSE)
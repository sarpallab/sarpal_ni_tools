# 2026-07-06 AndyP
# wrap_res_with_age.R
# wrapper function to get GAMs for MRSI

basedir <- '/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper/'
setwd(basedir)

source('/Users/andrew/Documents/GitHub/MRSI_gaba_glu/res_with_age.R')

df <- read_csv('13MP20200207_LCMv2fixidx_Raw.csv')
df <- df %>% separate_wider_delim(cols = ld8,delim="_",names=c("id","dateNumeric"),cols_remove=FALSE)
df$dateNumeric <- as.numeric(df$dateNumeric)


roi_list <- c('R Caudate','L Caudate','R Thalamus','L Thalamus')
met_list <- c('Cr','Asp.Cr','Cho.Cr','GABA.Cr','Glc.Cr','Gln.Cr','Glu.Cr','GPC.Cr','GSH.Cr','mI.Cr','GPC.Cho.Cr','NAA.NAAG.Cr','Glu.Gln.Cr','NAA.Cr','NAAG.Cr','Tau.Cr','Gpc.Cho.Cr')
met <- NULL
met_out <- NULL
M <- list()

for (iR in 1:length(roi_list)){
  df0 <- df %>% filter(label == roi_list[iR])
  for (iM in 1:length(met_list)){
    target_col0 <- paste0(met_list[iM],'.SD')
    target_col1 <- met_list[iM]
    df1 <- df %>% mutate(target_col1 = replace(met_list[iM], target_col0 > 20, NA))
    M0 <- res_with_age(MRSI_input = df1, met_name =met_list[iM],return_df=F,return_model = T, min_age = 18)
    M[[paste0("model_",roi_list[iR],'_',met_list[iM])]] <- M0
    target_col <- paste0(met_list[iM],'_gamadj')
    met0 <- res_with_age(MRSI_input = df, met_name =met_list[iM],return_df=T,return_model = F, min_age = 18) %>%
      select('id','dateNumeric','ld8','label',all_of(target_col))
    if (length(met)>0){
      met <- left_join(met,met0,by=c('id','dateNumeric','ld8','label'))
    } else {
      met <- met0
    }
  }
  print(roi_list[iR])
  met_out <- rbind(met_out,met)
}
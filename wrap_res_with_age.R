# 2026-07-06 AndyP
# wrap_res_with_age.R
# wrapper function to get GAMs for MRSI

library(tidyverse)

#basedir <- '/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper/'
basedir <- '/ix1/ginger/dsarpal/lab/reorg/projects/20260626-MRSI-Complete'
setwd(basedir)

#source('/Users/andrew/Documents/GitHub/MRSI_gaba_glu/res_with_age.R')
source('res_with_age.R')

df <- read_csv('13MP20200207_LCMv2fixidx_Raw.csv')
df <- df %>% separate_wider_delim(cols = ld8,delim="_",names=c("id","dateNumeric"),cols_remove=FALSE)
df$dateNumeric <- as.numeric(df$dateNumeric)
df <- df %>% group_by(id,visitnum,label) %>% slice(1) %>% ungroup()

roi_list <- c('R Caudate','L Caudate','R Thalamus','L Thalamus')
met_list <- c('Cr','Asp.Cr','Cho.Cr','GABA.Cr','Glc.Cr','Gln.Cr','Glu.Cr','GPC.Cr','GSH.Cr','mI.Cr','GPC.Cho.Cr','NAA.NAAG.Cr','Glu.Gln.Cr','NAA.Cr','NAAG.Cr','Tau.Cr')
met_out <- NULL
M <- list()

for (iR in 1:length(roi_list)){
  met <- NULL
  df0 <- df %>% filter(label == roi_list[iR])
  for (iM in 1:length(met_list)){
    target_col0 <- paste0(met_list[iM],'.SD')
    target_col1 <- met_list[iM]
    df1 <- df0 %>% mutate(target_col1 = replace(met_list[iM], target_col0 > 20, NA))
    df1 <- df0 %>% mutate(target_col1 = replace(met_list[iM], target_col1 > 5*sd(target_col1,na.rm=TRUE), NA))
    M0 <- res_with_age(MRSI_input = df1, met_name =met_list[iM],return_df=F,return_model = T, min_age = 18)
    model_name <- paste0("model_",roi_list[iR],'_',met_list[iM])
    M[[model_name]] <- M0
    #save(M0, file=paste0(model_name,'.Rdata'))
    target_col <- paste0(met_list[iM],'_gamadj')
    met0 <- res_with_age(MRSI_input = df1, met_name =met_list[iM],return_df=T,return_model = F, min_age = 18) %>%
      select('id','label','visitnum',all_of(target_col))
    if (length(met)>0){
      met <- full_join(met,met0,by=c('id','label','visitnum'))
    } else {
      met <- met0
    }
  }
  print(roi_list[iR])
  met_out <- rbind(met_out,met)
}

met_out1 <- inner_join(df,met_out,by=c('id','visitnum','label'))

save(M,file='20260706-gammodels.Rdata')
save(met_out1, file='20260706-gamadj-HC.Rdata')

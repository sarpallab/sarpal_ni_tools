# 2026-07-06 AndyP
# wrap_res_with_age.R
# wrapper function to get GAMs for MRSI

library(tidyverse)

basedir <- '/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper/'
#basedir <- '/ix1/ginger/dsarpal/lab/reorg/projects/20260626-MRSI-Complete'
setwd(basedir)

source('/Users/andrew/Documents/GitHub/MRSI_gaba_glu/res_with_age.R')
#source('res_with_age.R')

df <- read_csv('13MP20200207_LCMv2fixidx_Raw.csv')
df <- df %>% separate_wider_delim(cols = ld8,delim="_",names=c("id","dateNumeric"),cols_remove=FALSE)
df$dateNumeric <- as.numeric(df$dateNumeric)
df <- df %>% group_by(id,visitnum,label) %>% slice(1) %>% ungroup()

hc_mike <- read_excel('13MP20200207_LCMv2fixidx_Mike.xlsx') %>% separate_wider_delim(cols = RECID, delim = "_", names = c('id','date'))
#rm(met_out1)
#df <- df %>% filter(id %in% hc_mike$id)

roi_list <- c('R Caudate','L Caudate','R Thalamus','L Thalamus')
met_list <- c('Cr','Asp.Cr','Cho.Cr','GABA.Cr','Glc.Cr','Gln.Cr','Glu.Cr','GPC.Cr','GSH.Cr','mI.Cr','GPC.Cho.Cr','NAA.NAAG.Cr','Glu.Gln.Cr','NAA.Cr','NAAG.Cr','Tau.Cr')
met_out <- NULL
M <- list()

df <- df %>% mutate(Cr = case_when(Cr.SD > 20 | Cr < 0.01 | Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ Cr))
df <- df %>% mutate(Asp.Cr = case_when(Asp.SD > 20 | Asp.Cr < 0.01 | Asp.Cr > 5*sd(Asp.Cr,na.rm=T) ~ NA_real_, TRUE ~ Asp.Cr))
df <- df %>% mutate(Cho.Cr = case_when(Cho.SD > 20 | Cho.Cr < 0.01 | Cho.Cr > 5*sd(Cho.Cr,na.rm=T) ~ NA_real_, TRUE ~ Cho.Cr))
df <- df %>% mutate(GABA.Cr = case_when(GABA.SD > 20 | GABA.Cr < 0.01 | GABA.Cr > 5*sd(GABA.Cr,na.rm=T) ~ NA_real_, TRUE ~ GABA.Cr))
df <- df %>% mutate(Glc.Cr = case_when(Glc.SD > 20 | Glc.Cr < 0.01 | Glc.Cr > 5*sd(Glc.Cr,na.rm=T) ~ NA_real_, TRUE ~ Glc.Cr))
df <- df %>% mutate(Gln.Cr = case_when(Gln.SD > 20 | Gln.Cr < 0.01 | Gln.Cr > 5*sd(Gln.Cr,na.rm=T) ~ NA_real_, TRUE ~ Gln.Cr))
df <- df %>% mutate(Glu.Cr = case_when(Glu.SD > 20 | Glu.Cr < 0.01 | Glu.Cr > 5*sd(Glu.Cr,na.rm=T) ~ NA_real_, TRUE ~ Glu.Cr))
df <- df %>% mutate(GPC.Cr = case_when(GPC.SD > 20 | GPC.Cr < 0.01 | GPC.Cr > 5*sd(GPC.Cr,na.rm=T) ~ NA_real_, TRUE ~ GPC.Cr))
df <- df %>% mutate(GSH.Cr = case_when(GSH.SD > 20 | GSH.Cr < 0.01 | GSH.Cr > 5*sd(GSH.Cr,na.rm=T) ~ NA_real_, TRUE ~ GSH.Cr))
df <- df %>% mutate(mI.Cr = case_when(mI.SD > 20 | mI.Cr < 0.01 | mI.Cr > 5*sd(mI.Cr,na.rm=T) ~ NA_real_, TRUE ~ mI.Cr))
df <- df %>% mutate(GPC.Cho.Cr = case_when(GPC.Cho.SD > 20 | GPC.Cho.Cr < 0.01 | GPC.Cho.Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ GPC.Cho.Cr))
df <- df %>% mutate(NAA.NAAG.Cr = case_when(NAA.NAAG.SD > 20 | NAA.NAAG.Cr < 0.01 | NAA.NAAG.Cr > 5*sd(NAA.NAAG.Cr,na.rm=T) ~ NA_real_, TRUE ~ NAA.NAAG.Cr))
df <- df %>% mutate(Glu.Gln.Cr = case_when(Glu.Gln.SD > 20 | Glu.Gln.Cr < 0.01 | Glu.Gln.Cr > 5*sd(Glu.Gln.Cr,na.rm=T) ~ NA_real_, TRUE ~ Glu.Gln.Cr))
df <- df %>% mutate(NAA.Cr = case_when(NAA.SD > 20 | NAA.Cr < 0.01 | NAA.Cr > 5*sd(NAA.Cr,na.rm=T) ~ NA_real_, TRUE ~ NAA.Cr))
df <- df %>% mutate(NAAG.Cr = case_when(NAAG.SD > 20 | NAAG.Cr < 0.01 | NAAG.Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ NAAG.Cr))
df <- df %>% mutate(Tau.Cr = case_when(Tau.SD > 20 | Tau.Cr < 0.01 | Tau.Cr > 5*sd(Tau.Cr,na.rm=T) ~ NA_real_, TRUE ~ Tau.Cr))


for (iR in 1:length(roi_list)){
  met <- NULL
  df0 <- df %>% filter(label == roi_list[iR])
  for (iM in 1:length(met_list)){
    
    # if (iR == 3){
    #   if (met_list[iM] == 'NAAG.Cr'){
    #     browser()
    #   }
    # }
    
    target_col0 <- paste0(met_list[iM],'.SD')
    target_col1 <- met_list[iM]
    M0 <- res_with_age(MRSI_input = df0, met_name =met_list[iM],return_df=F,return_model = T, min_age = 18)
    model_name <- paste0("model_",roi_list[iR],'_',met_list[iM])
    M[[model_name]] <- M0
    #save(M0, file=paste0(model_name,'.Rdata'))
    target_col <- paste0(met_list[iM],'_gamadj')
    met0 <- res_with_age(MRSI_input = df0, met_name =met_list[iM],return_df=T,return_model = F, min_age = 18) %>%
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

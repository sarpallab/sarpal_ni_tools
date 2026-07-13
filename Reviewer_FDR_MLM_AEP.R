# 2026-06-03 AndyP
# address reviewer comments:

# 2026-06-17 AndyP
# added contralateral thalamus

library(tidyverse)
library(readxl)
library(lmerTest)
library(broom.mixed)
library(cocor)
library(emmeans)
library(DescTools)
library(ggrepel)
library(ggsignif)

reload_new_20260706 = T# 2026-07-06 reload after redoing gamadj models AndyP
nan_out_Crgamadj = F # 2026-07-13 AndyP.  Subsetting reload_new_20260706, nan-ning out Cr_gamadj where SD > 20, values are above the floor, and below 5 * mean values.
# This results in Ca Glu N.S. on right side, might still be useful sensitivity analysis??
reload_new = F # 2026-07-02 post redoing GM
# same as reload, but loads in non gam-adjusted metabolites, can merge with reload df manually
# also does group analyses with non-gam-adjusted metabolites
# old, part of redoing GM for revision using uniform method
reload = F
reload_check_GM = F
longitudinal = F
group = F
hilowdoi = F
eibalance = F
group_r_nr = F
# will reload just Sarpal / CZ data
clinical = F
handedness_group = F
Figure_2 = F
Creatine_Check  = T
loglink_GLM = F
Figure_3 = F
Figure_4 = F
corr_heatmap = F
#The manuscript states that FDR correction was performed by accounting for metabolites within each ROI. 
# However, the statistical inference and biological interpretation are made across three ROIs. 
# Please state exactly which p-values were included in each FDR correction family. 
# In addition, please provide a more conservative sensitivity analysis in which FDR correction 
# includes all ROI-by-metabolite tests within each type of analysis, and report whether the main findings remain significant.

# macbook
 #basedir <- ('/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper')
# macmini
basedir <- '/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper'
setwd(basedir)

#df <- read_csv('20260708-final-dataset-MRSI-2.csv')

if (reload_new_20260706 == T){
  load('20260706-gamadj-HC.Rdata') # loads met_out1
  hc <- met_out1
  hc_mike <- read_excel('13MP20200207_LCMv2fixidx_Mike.xlsx') %>% separate_wider_delim(cols = RECID, delim = "_", names = c('id','date'))
  rm(met_out1)
  hc <- hc %>% filter(id %in% hc_mike$id)
  
  
  gc()
  load('20260706-SSD-gamadj.Rdata') # loads sz_met_out
  ssd <- sz_met_out %>% group_by(RECID,timepoint,region,metabolite) %>% slice(1) %>% ungroup()
  rm(sz_met_out)
  
  hc <- hc %>% select(id,visitnum,sex,age,label,GMrat,all_of(contains('_gamadj')))
  hc <- hc %>% rename(timepoint = visitnum)
  
  ssd <- ssd %>% pivot_wider(id_cols = c(RECID,timepoint,region,age), names_from = metabolite, values_from = value)
  
  # ssd did not have sex
  sd <- read_excel('7T.DARES.baselines.xlsx', sheet = 3) %>% group_by(RECID) %>% slice(1) %>% ungroup() %>% select(RECID,sex,"DUP (months)" )
  sd <- sd %>% mutate(sex = as.factor(case_when(sex == 1 ~ 'M',
                                                sex == 2 ~ 'F')))
  
  ssd <- inner_join(ssd,sd,by='RECID')
  
  # get GMrat
  # note that Mike used 
  
  ssdgm <- read_excel('sarpal_mrsi_original_07062026.xlsx') %>% select(RECID,timepoint,region,GMrat) %>% group_by(RECID,timepoint,region) %>% slice(1) %>% ungroup()
  ssdgm <- ssdgm %>% mutate(region = case_when(region == 'R Caudate' ~ 'R Caudate',
                                                 region == 'right caudate' ~ 'R Caudate',
                                                 region == 'L Caudate' ~ 'L Caudate',
                                                 region == 'left caudate' ~ 'L Caudate',
                                                 region == 'R Thalamus' ~ 'R Thalamus',
                                                 region == 'right thalamus' ~ 'R Thalamus',
                                                 region == 'L Thalamus' ~ 'L Thalamus',
                                                 region == 'left thalamus' ~ 'L Thalamus'))
  
  ssd <- inner_join(ssd,ssdgm,by=c('RECID','timepoint','region'))
  
  ssd <- ssd %>% rename(id = RECID)
  hc <- hc %>% rename(region = label)
  hc <- hc %>% mutate("DUP (months)"  = NA)
  # to check, add SD values, check that results hold when covarying for CRLB, can also check FD
  
  ssd <- ssd %>% mutate(group_level = 'SZ')
  hc <- hc %>% mutate(group_level = 'HC')
  
  common_col <- intersect(colnames(hc),colnames(ssd))
  hc <- hc %>% select(common_col)
  ssd <- ssd %>% select(common_col)
  
  df <- rbind(ssd,hc)
  
  ############################
  #add in remitter status from outside spreadsheet
  supplemental_data2 <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_data2 <- supplemental_data2 %>% filter(Timepoint != 0)
  supplemental_data2 <- supplemental_data2 %>%
    mutate(Remitter_Status = ifelse(CONCDIS4 <= 3 & HALL12 <= 3 & UTC15 <= 3, "R", "NR"), timepoint = Timepoint)
  supplemental_data2 <- supplemental_data2 %>% select(RECID,POSSX,Remitter_Status, timepoint)
  supplemental_dataBL <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_dataBL <- supplemental_dataBL %>% filter(Timepoint == 0) %>% mutate(timepoint = Timepoint, Remitter_Status = NA)
  supplemental_dataBL <- supplemental_dataBL %>% select(RECID, POSSX,timepoint,Remitter_Status)
  supplemental_data2 <- rbind(supplemental_data2, supplemental_dataBL)
  supplemental_data2 <- supplemental_data2 %>% group_by(RECID) %>% tidyr::fill(Remitter_Status, .direction = "downup") %>% ungroup() %>% mutate(timepoint = case_when(timepoint == 0 ~ 1,
                                                                                                                                                                      timepoint > 0 ~ 2))
  # sarpal_data_BLFU <- sarpal_data_BLFU %>% 
  #   mutate(Remitter_Status = supplemental_data2$Remitter_Status[match(RECID,supplemental_data2$RECID)])
  supplemental_data2$RECID <- as.character(supplemental_data2$RECID)
  supplemental_data2 <- supplemental_data2 %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                                                            timepoint == 2 ~ 'FU')) %>% rename(id = RECID)
  
  df <- df %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                            timepoint == 2 ~ 'FU',
                                            timepoint == 3 ~ 'extra'))
  df <- full_join(df,supplemental_data2,by=c('id','timepoint'))
  df$doi_m <- df$`DUP (months)`
  df <- df %>% dplyr::select(id,doi_m,Cr_gamadj, region,group_level,timepoint,age,sex,POSSX,Remitter_Status,GMrat,GPC.Cr_gamadj,Glc.Cr_gamadj,Glu.Cr_gamadj,GPC.Cho.Cr_gamadj,GABA.Cr_gamadj,NAA.Cr_gamadj,mI.Cr_gamadj,Gln.Cr_gamadj,NAAG.Cr_gamadj,Glu.Gln.Cr_gamadj)
  df <- df %>% mutate(roi = case_when(region == 'R Caudate' ~ 'R Caudate',
                                      region == 'right caudate' ~ 'R Caudate',
                                      region == 'L Caudate' ~ 'L Caudate',
                                      region == 'left caudate' ~ 'L Caudate',
                                      region == 'R Thalamus' ~ 'R Thalamus',
                                      region == 'right thalamus' ~ 'R Thalamus',
                                      region == 'L Thalamus' ~ 'L Thalamus',
                                      region == 'left thalamus' ~ 'L Thalamus'))
  df <- df %>% dplyr::mutate(hemi = str_sub(roi, 1,1),roi = str_sub(roi, 3, -1))
  
  
  # this sheet has NaN that were manually removed by Mike.  These NaN's need to be re-inserted into the GAM model outputs, which seems to have interpolated over the NaNs
  # check that values are identical for both BL and BL/FU groups
  sz_mets <- c('Cre','GABA/Cre', 'Glu/Cre', 'Gln/Cre','Glc/Cre','mI/Cre','NAA/Cre', 'Glu.Gln/Cre', 'NAAG/Cre', 'GPC/Cre', 'GPC.Cho/Cre')
  ndistinct <- read_excel('sarpal_mrsi_original_07062026.xlsx') %>% select(RECID,timepoint,region,GMrat) %>% group_by(RECID,timepoint,region) %>% 
    summarize(nD1 = n_distinct('GABA/Cre'),nD2 = n_distinct('Glu/Cre'),nD3 = n_distinct('Gln/Cre'),nD4 = n_distinct('Glc/Cre'),nD5 = n_distinct('mI/Cre'),nD6 = n_distinct('NAA/Cre'),nD7 = n_distinct('Glu.Gln/Cre'),nD8 = n_distinct('NAAG/Cre'),nD9 = n_distinct('GPC/Cre'),nD10 = n_distinct('GPC.Cho/Cre')) %>% ungroup()
  # all nD values are == 1, so we are good.
  ssdna <- read_excel('sarpal_mrsi_original_07062026.xlsx') %>% select(RECID,timepoint,region,all_of(sz_mets)) %>% group_by(RECID,timepoint,region) %>% slice(1) %>% ungroup()
  ssdna <- ssdna %>% mutate(region = case_when(region == 'R Caudate' ~ 'R Caudate',
                                               region == 'right caudate' ~ 'R Caudate',
                                               region == 'L Caudate' ~ 'L Caudate',
                                               region == 'left caudate' ~ 'L Caudate',
                                               region == 'R Thalamus' ~ 'R Thalamus',
                                               region == 'right thalamus' ~ 'R Thalamus',
                                               region == 'L Thalamus' ~ 'L Thalamus',
                                               region == 'left thalamus' ~ 'L Thalamus'))
  
  ssdna <- ssdna %>% arrange(RECID,timepoint,region)
  
  na_indices_GABA <- which(is.na(as.numeric(ssdna$`GABA/Cre`)), arr.ind = TRUE)
  na_indices_Glu <- which(is.na(as.numeric(ssdna$`Glu/Cre`)), arr.ind = TRUE)  
  na_indices_Gln <- which(is.na(as.numeric(ssdna$`Gln/Cre`)), arr.ind = TRUE)
  na_indices_Glc <- which(is.na(as.numeric(ssdna$`Glc/Cre`)), arr.ind = TRUE)
  na_indices_mI <- which(is.na(as.numeric(ssdna$`mI/Cre`)), arr.ind = TRUE)
  na_indices_NAA <- which(is.na(as.numeric(ssdna$`NAA/Cre`)), arr.ind = TRUE)
  na_indices_Glu.Gln <- which(is.na(as.numeric(ssdna$`Glu.Gln/Cre`)), arr.ind = TRUE)
  na_indices_NAAG <- which(is.na(as.numeric(ssdna$`NAAG/Cre`)), arr.ind = TRUE)
  na_indices_GPC <- which(is.na(as.numeric(ssdna$`GPC/Cre`)), arr.ind = TRUE)
  na_indices_GPC.Cho <- which(is.na(as.numeric(ssdna$`GPC.Cho/Cre`)), arr.ind = TRUE)

  dfssd <- df %>% filter(group_level == 'SZ') %>% arrange(id,timepoint,roi)
  dfhc <- df %>% filter(group_level == 'HC')
  
  dfssd$GABA.Cr_gamadj[na_indices_GABA] = NA
  dfssd$Glu.Cr_gamadj[na_indices_Glu] = NA
  dfssd$Gln.Cr_gamadj[na_indices_Gln] = NA
  dfssd$Glc.Cr_gamadj[na_indices_Glc] = NA
  dfssd$mI.Cr_gamadj[na_indices_mI] = NA
  dfssd$NAA.Cr_gamadj[na_indices_NAA] = NA
  dfssd$Glu.Gln.Cr_gamadj[na_indices_Glu.Gln] = NA
  dfssd$GPC.Cr_gamadj[na_indices_GPC] = NA
  dfssd$GPC.Cho.Cr_gamadj[na_indices_GPC.Cho] = NA
  # testing 2/3 rule (2/3 of data is NA => exclude metabolite)
  
  df <- rbind(dfssd,dfhc)
  
  test_in <- df %>% 
    group_by(roi, group_level, timepoint,hemi) %>% 
    summarize(
      N = n_distinct(id),
      # Use lowercase n() and make sure 'gamadj' matches your column names exactly
      across(contains('gamadj'), ~ sum(is.na(.)) / n()), 
      .groups = "drop" # Replaces the need for a separate ungroup()
    )
  
  df <- df %>% 
    group_by(roi, group_level, timepoint,hemi) %>% 
    mutate(across(
      contains("gamadj"), 
      ~ case_when(
        mean(is.na(.)) >= 1/3 ~ NA,
        TRUE ~ . 
      )
    )) %>% 
    ungroup()
  
  test_out <- df %>% 
    group_by(roi, group_level, timepoint,hemi) %>% 
    summarize(
      N = n_distinct(id),
      # Use lowercase n() and make sure 'gamadj' matches your column names exactly
      across(contains('gamadj'), ~ sum(is.na(.)) / n()), 
      .groups = "drop" # Replaces the need for a separate ungroup()
    )
  
  df <- df %>% dplyr::rename(id = 'id',GPC = 'GPC.Cr_gamadj',Glu = 'Glu.Cr_gamadj',GPC.Cho = 'GPC.Cho.Cr_gamadj',Glc = 'Glc.Cr_gamadj',GABA = 'GABA.Cr_gamadj',NAA = 'NAA.Cr_gamadj',mI = 'mI.Cr_gamadj',Gln = 'Gln.Cr_gamadj',NAAG = 'NAAG.Cr_gamadj',Glu.Gln = 'Glu.Gln.Cr_gamadj')
  
  df <- df %>% filter(age >= 18)
  
  df <- df %>% group_by(id,roi,hemi) %>% mutate(nsess = 1:n()) %>%
    ungroup()
  
  df$timepoint <- as.numeric(df$timepoint)
  
  
  # 11725, 11810 use extra
  df <- df %>% group_by(id) %>% arrange(nsess) %>% 
    mutate(timepoint1 = case_when(nsess == 1 ~ 'BL',
                                  nsess == 2 ~ 'FU',
                                  nsess >= 3 ~ 'discard')
    )
  df <- df %>% filter(timepoint1 != 'discard') %>%
    select(!timepoint) %>% rename(timepoint = timepoint1)
  
  # this should not reduce the dataframe any further
  #df <- df %>% group_by(id,timepoint,roi,hemi) %>% slice(1) %>% ungroup()
  
  df$timepoint <- relevel(factor(df$timepoint), ref = "BL")
  df$sex <- relevel(factor(df$sex),ref = "F")
  df <- df %>% mutate(age_sc = scale(age))
  
  df <- df %>% mutate(group = case_when(group_level == 'HC' ~ 'HC',
                                        Remitter_Status == 'NR' ~ 'NR',
                                        Remitter_Status == 'R' ~ 'R'))

  df$group <- relevel(factor(df$group), ref = "HC")
  df <- df %>%
    mutate(
      # Combine Group and Time into 3 valid conditions
      condition = case_when(
        group_level == "HC" ~ "HC_BL",
        Remitter_Status == "NR" & timepoint == "BL"  ~ "NR_BL",
        Remitter_Status == "NR" & timepoint == "FU"  ~ "NR_FU",
        Remitter_Status == "R" & timepoint == "BL"  ~ "R_BL",
        Remitter_Status == "R" & timepoint == "FU"  ~ "R_FU",      
      ),
      condition = factor(condition, levels = c("HC_BL", "NR_BL", "R_BL", "NR_FU","R_FU"))
    )
  
  df$condition <- relevel(df$condition, ref = "HC_BL")
  
  df$group_level <- relevel(as.factor(df$group_level),'HC')
  df0 <- df %>% group_by(id) %>% slice(1) %>% ungroup()
  mDOI <- median(df0$doi_m[df0$group_level == 'SZ'],na.rm=TRUE)
  df <- df %>% mutate(median_split_doi = case_when(doi_m >= mDOI ~ 'high',
                                                   doi_m < mDOI ~ 'low')
  )
  df$median_split_doi[df$id == "2695"] = 'low'
  
  # df <- df %>% mutate(Glu = case_when(Glu < 5*mean(Glu,na.rm=T) ~ Glu,
  #                                     Glu >= 5*mean(Glu,na.rm=T) ~ NA))
  # df <- df %>% mutate(NAAG = case_when(NAAG <5*mean(NAAG,na.rm=T) ~ NAAG,
  #                                      NAAG >= 5*mean(NAAG,na.rm=T) ~ NA))
  # df <- df %>% mutate(NAA = case_when(NAA < 5*mean(NAA,na.rm=T) ~ NAA,
  #                                      NAA >= 5*mean(NAA,na.rm=T) ~ NA))
  # df <- df %>% mutate(mI = case_when(mI < 5*mean(mI,na.rm=T) ~ mI,
  #                                      mI >= 5*mean(mI,na.rm=T) ~ NA))
  # df <- df %>% mutate(Glu.Gln = case_when(Glu.Gln < 5*mean(Glu.Gln,na.rm=T) ~ Glu.Gln,
  #                                      Glu.Gln >= 5*mean(Glu.Gln,na.rm=T) ~ NA))
  
  
  dg <- read_csv('13MP20200207_LCMv2fixidx_Raw.csv')
  dg <- dg %>% separate_wider_delim(cols = ld8,delim="_",names=c("id","dateNumeric"),cols_remove=FALSE)
  dg$dateNumeric <- as.numeric(dg$dateNumeric)
  dg <- dg %>% group_by(id,visitnum,label) %>% slice(1) %>% ungroup()
  dg <- dg %>% mutate(Cr = case_when(Cr.SD > 20 | Cr < 0.01 | Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ Cr))
  dg <- dg %>% filter(age >= 18) %>%
    filter(label %in% c('L Thalamus','R Thalamus','L Caudate','R Caudate')) %>% 
    group_by(id,label) %>% mutate(nsess = 1:n()) %>%
    filter(nsess == 1) %>% mutate(timepoint = 'BL')
  
  
  # Match ROI names, create date value
  dg <- dg %>% 
    separate(label, c('hemi','roi'), ' ', convert=TRUE) %>%
    mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
    mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
    mutate(hemi = ifelse(hemi=='L', 'L', ifelse(hemi=='R', 'R', NA))) %>% # Fixed here
    mutate(region = paste0(hemi,' ', roi))
  dg <- dg %>% select(id,roi,hemi,timepoint,Cr, Cr.SD)

  
if (nan_out_Crgamadj==T){  
    
  hc_mike <- read_excel('13MP20200207_LCMv2fixidx_Mike.xlsx') %>% separate_wider_delim(cols = RECID, delim = "_", names = c('id','date'))
  #rm(met_out1)
  #Load SZ MRSI data
  szmet_orig <- readxl::read_xlsx('sarpal_mrsi_original_07062026.xlsx', sheet = 1) %>%
    rename(scan_date = Scan_date) %>% mutate(source = 'orig')
  #szmet_orig <- szmet_orig %>% group_by(RECID) %>% separate_wider_delim(cols = region, delim = " ",names = c('hemisphere','roi'), cols_remove = TRUE,too_few = "align_start") %>% ungroup()
  szmet_orig <- szmet_orig %>% rename(roi = 'region')
  # Load SZ MRSI data
  szmet_orig <- readxl::read_xlsx('sarpal_mrsi_original_07062026.xlsx', sheet = 1) %>%
    rename(scan_date = Scan_date) %>% mutate(source = 'orig')
  #szmet_orig <- szmet_orig %>% group_by(RECID) %>% separate_wider_delim(cols = region, delim = " ",names = c('hemisphere','roi'), cols_remove = TRUE,too_few = "align_start") %>% ungroup()
  szmet_orig <- szmet_orig %>% rename(roi = 'region') %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                                                                       timepoint == 2 ~ 'FU'))


  # Match ROI names, create date value
  szmet_new <- szmet_orig %>%
    separate(roi, c('hemi','roi'), ' ', convert=TRUE) %>%
    mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
    mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
    mutate(dateNumeric = as.numeric(as.POSIXct(scan_date, format="%Y-%m-%d")),
           hemi = ifelse(hemi=='left', 'L', ifelse(hemi=='right', 'R', NA))) %>%
    mutate(region = paste0(hemi,' ', roi))

  # Match ROI names, create date value
  szmet_new <- szmet_orig %>%
    separate(roi, c('hemi','roi'), ' ', convert=TRUE) %>%
    mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
    mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
    mutate(dateNumeric = as.numeric(as.POSIXct(scan_date, format="%Y-%m-%d")),
           hemi = ifelse(hemi=='left', 'L', ifelse(hemi=='right', 'R', NA))) %>%
    mutate(region = paste0(hemi,' ', roi))

  szmet_new <- szmet_new %>% select(RECID,timepoint,roi,hemi,Cre,`Cre %SD`)
  szmet_new <- szmet_new %>% rename(id = RECID, Cr = Cre, Cr.SD = `Cre %SD`)
  szmet_new$id <- as.character(szmet_new$id)
  szmet_new <- szmet_new %>% mutate(Cr = case_when(Cr.SD > 20 | Cr < 0.01 | Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ Cr))

  df <- df %>% filter((group_level == 'HC' & timepoint == 'BL') | group_level == 'SZ')


  dg <- rbind(dg,szmet_new)
  df <- left_join(df,dg,by=c('id','timepoint','roi','hemi'))

  df <- df %>% mutate(Cr_gamadj = case_when(Cr.SD > 20 | Cr < 0.01 | Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ Cr))

  na_indices_Cr <- which(is.na(as.numeric(df$Cr_gamadj)), arr.ind = TRUE)

  df$GABA[na_indices_Cr] = NA
  df$Glu[na_indices_Cr] = NA
  df$Gln[na_indices_Cr] = NA
  df$Glc[na_indices_Cr] = NA
  df$mI[na_indices_Cr] = NA
  df$NAA[na_indices_Cr] = NA
  df$Glu.Gln[na_indices_Cr] = NA
  df$GPC[na_indices_Cr] = NA
  df$GPC.Cho[na_indices_Cr] = NA
}
  
  df <- df %>% select(!Glc) # No Glc data in SSD
  df0 <- df %>% group_by(group_level,id) %>% slice(1) %>% ungroup() %>% group_by(group_level) %>% summarize(mA = mean(age, na.rm=TRUE), sA = sd(age,na.rm=TRUE), N= n()) %>% ungroup()
  dfbprs <- df %>% filter(group_level == 'SZ') %>% group_by(id,timepoint) %>% slice(1) %>% ungroup() %>% group_by(timepoint) %>% summarize(mbprs = mean(POSSX,na.rm=TRUE), sbprs = sd(POSSX,na.rm=TRUE)) %>% ungroup()
  
  df2 <- df %>% group_by(group_level,id) %>% slice(1) %>% ungroup()
  df1 <- df2 %>% select(group_level,sex) %>% mutate(sex = case_when(sex == 'M' ~ 1, sex == 'F' ~ 2), group_level = case_when(group_level == 'HC' ~ 1, group_level == 'SZ' ~ 2))
  df1 <- as.matrix(df1)
  counts_table <- table(df1[, "group_level"], df1[, "sex"])
  st <- chisq.test(counts_table)
  
}


if (reload_new == T){
  hc <- read_csv(paste0(basedir,'/HC_20260703/','metabolites_gamadj_HC_20260703.csv')) %>% mutate(group = 'HC', timepoint = 'BL',`DUP (months)` = NA)
  ssd <- read_csv(paste0(basedir,'/SSD/','metabolites_gamadj_long_07012026.csv')) %>% mutate(group = 'NA',
                                                                                             timepoint = case_when(timepoint == 1 ~ 'BL',
                                                                                                                   timepoint == 2 ~ 'FU'))
  # ssd did not have sex
  sd <- read_excel('7T.DARES.baselines.xlsx', sheet = 3) %>% group_by(RECID) %>% slice(1) %>% ungroup() %>% select(RECID,sex,"DUP (months)" )
  sd <- sd %>% mutate(sex = as.factor(case_when(sex == 1 ~ 'M',
                                                sex == 2 ~ 'F')))
  
  ssd <- inner_join(ssd,sd,by='RECID')
  
  common_cols <- intersect(colnames(hc),colnames(ssd))
  hc <- hc %>% select(all_of(common_cols))
  ssd <- ssd %>% select(all_of(common_cols))
  df <- rbind(hc,ssd) 
  
  mike <- read_xlsx('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper/HC_20260703/13MP20200207_LCMv2fixidx_Mike.xlsx')
  mike <- mike %>% select("RECID","GMrat","Glu.Cr","age") %>% rename(GMrat_mike = GMrat, Glu.Cr_mike = Glu.Cr) %>% filter(age >= 18) %>% group_by(RECID) %>% slice(1) %>% ungroup() %>% select(!age)
  
  df <- left_join(df,mike,by=c('RECID'))
  
  
  ############################
  #add in remitter status from outside spreadsheet
  supplemental_data2 <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_data2 <- supplemental_data2 %>% filter(Timepoint != 0)
  supplemental_data2 <- supplemental_data2 %>%
    mutate(Remitter_Status = ifelse(CONCDIS4 <= 3 & HALL12 <= 3 & UTC15 <= 3, "R", "NR"), timepoint = Timepoint)
  supplemental_data2 <- supplemental_data2 %>% select(RECID,POSSX,Remitter_Status, timepoint)
  supplemental_dataBL <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_dataBL <- supplemental_dataBL %>% filter(Timepoint == 0) %>% mutate(timepoint = Timepoint, Remitter_Status = NA)
  supplemental_dataBL <- supplemental_dataBL %>% select(RECID, POSSX,timepoint,Remitter_Status)
  supplemental_data2 <- rbind(supplemental_data2, supplemental_dataBL)
  supplemental_data2 <- supplemental_data2 %>% group_by(RECID) %>% tidyr::fill(Remitter_Status, .direction = "downup") %>% ungroup() %>% mutate(timepoint = case_when(timepoint == 0 ~ 1,
                                                                                                                                                                      timepoint > 0 ~ 2))
  # sarpal_data_BLFU <- sarpal_data_BLFU %>% 
  #   mutate(Remitter_Status = supplemental_data2$Remitter_Status[match(RECID,supplemental_data2$RECID)])
  supplemental_data2$RECID <- as.character(supplemental_data2$RECID)
  supplemental_data2 <- supplemental_data2 %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                                                            timepoint == 2 ~ 'FU'))
  df <- left_join(df,supplemental_data2,by=c('RECID','timepoint'))
  df$doi_m <- df$`DUP (months)`
  df <- df %>% dplyr::select(RECID,doi_m,region,timepoint,age,sex,group,Scan_date,POSSX,Remitter_Status,GMrat,GMrat_mike,Glu.Cr_mike,GPC.Cr_gamadj,Glc.Cr_gamadj,Glu.Cr_gamadj,GPC.Cho.Cr_gamadj,GABA.Cr_gamadj,NAA.Cr_gamadj,mI.Cr_gamadj,Gln.Cr_gamadj,NAAG.Cr_gamadj,Glu.Gln.Cr_gamadj,Glu.Cr)
  df <- df %>% mutate(roi = case_when(region == 'R Caudate' ~ 'R Caudate',
                                      region == 'right caudate' ~ 'R Caudate',
                                      region == 'L Caudate' ~ 'L Caudate',
                                      region == 'left caudate' ~ 'L Caudate',
                                      region == 'R Thalamus' ~ 'R Thalamus',
                                      region == 'right thalamus' ~ 'R Thalamus',
                                      region == 'L Thalamus' ~ 'L Thalamus',
                                      region == 'left thalamus' ~ 'L Thalamus'))
  df <- df %>% dplyr::mutate(hemi = str_sub(roi, 1,1),roi = str_sub(roi, 3, -1))
  df <- df %>% 
    group_by(roi, group, timepoint) %>% 
    mutate(across(
      contains("gamadj"), 
      ~ case_when(
        mean(is.na(.)) >= 2/3 ~ NA,
        TRUE ~ . 
      )
    )) %>% 
    ungroup()
  df <- df %>% dplyr::rename(id = 'RECID',GPC = 'GPC.Cr_gamadj',Glu = 'Glu.Cr_gamadj',GPC.Cho = 'GPC.Cho.Cr_gamadj',Glc = 'Glc.Cr_gamadj',GABA = 'GABA.Cr_gamadj',NAA = 'NAA.Cr_gamadj',mI = 'mI.Cr_gamadj',Gln = 'Gln.Cr_gamadj',NAAG = 'NAAG.Cr_gamadj',Glu.Gln = 'Glu.Gln.Cr_gamadj')
  df$timepoint <- relevel(factor(df$timepoint), ref = "BL")
  df$sex <- relevel(factor(df$sex),ref = "F")
  df <- df %>% mutate(age_sc = scale(age))
  
  df <- df %>% mutate(group = case_when(group == 'HC' ~ 'HC',
                                        Remitter_Status == 'NR' ~ 'NR',
                                        Remitter_Status == 'R' ~ 'R')) %>%
    select(!Remitter_Status)
  
  df <- df %>% mutate(group_level = case_when(group == 'HC' ~ 'HC',
                                              group == 'NR' ~ 'SZ',
                                              group == 'R' ~ 'SZ',
                                              is.na(group) ~ 'SZ')
  )
  df$group <- relevel(factor(df$group), ref = "HC")
  df <- df %>%
    mutate(
      # Combine Group and Time into 3 valid conditions
      condition = case_when(
        group == "HC" ~ "HC_BL",
        group == "NR" & timepoint == "BL"  ~ "NR_BL",
        group == "NR" & timepoint == "FU"  ~ "NR_FU",
        group == "R" & timepoint == "BL"  ~ "R_BL",
        group == "R" & timepoint == "FU"  ~ "R_FU",      
      ),
      condition = factor(condition, levels = c("HC_BL", "NR_BL", "R_BL", "NR_FU","R_FU"))
    )
  
  df$condition <- relevel(df$condition, ref = "HC_BL")
  
  df$group_level <- relevel(as.factor(df$group_level),'HC')
  df0 <- df %>% group_by(id) %>% slice(1) %>% ungroup()
  mDOI <- median(df0$doi_m[df0$group_level == 'SZ'],na.rm=TRUE)
  df <- df %>% mutate(median_split_doi = case_when(doi_m >= mDOI ~ 'high',
                                                   doi_m < mDOI ~ 'low')
  )
  df$median_split_doi[df$id == "2695"] = 'low'
  
  df <- df %>% mutate(Glu = case_when(Glu < 7.5 ~ Glu,
                                      Glu >= 7.5 ~ NA))
  
  df <- df %>% filter(age >= 18)
  
  df0 <- df %>% group_by(group_level,id) %>% slice(1) %>% ungroup() %>% group_by(group_level) %>% summarize(mA = mean(age, na.rm=TRUE), sA = sd(age,na.rm=TRUE), N= n()) %>% ungroup()
  dfbprs <- df %>% filter(group_level == 'SZ') %>% group_by(id,timepoint) %>% slice(1) %>% ungroup() %>% group_by(timepoint) %>% summarize(mbprs = mean(POSSX,na.rm=TRUE), sbprs = sd(POSSX,na.rm=TRUE)) %>% ungroup()
  df <- df %>% mutate(Glu.Cr_mike = case_when(group_level == 'HC' ~ Glu.Cr_mike,
                                              group_level == 'SZ' ~ Glu),
                      GMrat_mike = case_when(group_level == 'HC' ~ GMrat_mike,
                                             group_level == 'SZ' ~ GMrat))
    
}

if (reload_check_GM==T){
  
  #read in ALL corrected data (including duplicates)
  #isolate patients with baseline and follow-up only (ie only patients with known remitter status)
  sarpal_data_all <- read_csv('sarpal_mrsi_adj_all_20250310.csv') %>%
    filter(source %in% c('BL.FU','DARES.baselines')) %>% group_by(RECID,timepoint,hemi,roi) %>% slice(1) %>% ungroup()
  
  #add in sex from outside spreadsheet
  # 2026-06-02 AndyP changed to sheet 3 was reading incorrect sheet
  sarpal_data_supplemental <- read_excel('7T.DARES.baselines.xlsx', sheet = 3) %>% group_by(RECID) %>% slice(1) %>% ungroup()
  
  df <- inner_join(sarpal_data_all,sarpal_data_supplemental,by='RECID')
  
  df <- df %>% mutate(sex = as.factor(case_when(sex == 1 ~ 'M',
                                                sex == 2 ~ 'F')))
  
  
  ############################
  #add in remitter status from outside spreadsheet
  supplemental_data2 <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_data2 <- supplemental_data2 %>% filter(Timepoint != 0)
  supplemental_data2 <- supplemental_data2 %>%
    mutate(Remitter_Status = ifelse(CONCDIS4 <= 3 & HALL12 <= 3 & UTC15 <= 3, "R", "NR"), timepoint = Timepoint)
  supplemental_data2 <- supplemental_data2 %>% select(RECID,POSSX,Remitter_Status, timepoint)
  supplemental_dataBL <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_dataBL <- supplemental_dataBL %>% filter(Timepoint == 0) %>% mutate(timepoint = Timepoint, Remitter_Status = NA)
  supplemental_dataBL <- supplemental_dataBL %>% select(RECID, POSSX,timepoint,Remitter_Status)
  supplemental_data2 <- rbind(supplemental_data2, supplemental_dataBL)
  supplemental_data2 <- supplemental_data2 %>% group_by(RECID) %>% tidyr::fill(Remitter_Status, .direction = "downup") %>% ungroup() %>% mutate(timepoint = case_when(timepoint == 0 ~ 1,
                                                                                                                                                                      timepoint > 0 ~ 2))
  # sarpal_data_BLFU <- sarpal_data_BLFU %>% 
  #   mutate(Remitter_Status = supplemental_data2$Remitter_Status[match(RECID,supplemental_data2$RECID)])
  
  df <- left_join(df,supplemental_data2,by=c('RECID','timepoint'))
  df <- df %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                            timepoint == 2 ~ 'FU'))
  
  df$doi_m <- df$`DUP (months)`
  # note Cr_gamadj is Metabolite / Creatine NOT raw Creatine metabolite.
  # The adjustment is based on the length of the scan and was done for luna lab and Sarpal lab.
  load('luna_data.Rdata')
  #get rid of unimportant columns for our immediate analysis to allow for wide pivot
  luna_data <- luna_data[, c("met", "ld8", "region", "age", "GMrat","Cr_gamadj")]
  #pivot wide for easy comparison with other spreadsheet
  luna_data <- luna_data %>%
    pivot_wider(names_from = met, values_from = Cr_gamadj)
  luna_data <- luna_data %>%
    separate(ld8, into = c("lunaid", "date"), sep = "_")
  luna_data$lunaid <- as.numeric(luna_data$lunaid)
  luna_data <- luna_data %>% filter(age >= 18)
  luna_data <- luna_data %>% dplyr::rename(GPC0 = 'GPC',Glu0 = 'Glu',GPC.Cho0 = 'GPC.Cho',Glc0 = 'Glc',GABA0 = 'GABA',NAA0 = 'NAA',mI0 = 'mI',Gln0 = 'Gln',NAAG0 = 'NAAG',Glu.Gln0 = 'Glu.Gln')
  
  #add in sex from outside spreadsheet
  luna_data_supplemental <- read_csv('luna7t_sifpc_subset.csv')
  luna_data_supplemental <- luna_data_supplemental %>% dplyr::select(lunaid,sex) %>% group_by(lunaid) %>% slice(1)
  
  luna_data <- inner_join(luna_data,luna_data_supplemental,by=c("lunaid"))
  
  luna_data <- luna_data %>%
    mutate(date = as.numeric(date)) %>%
    arrange(region, lunaid, date) %>%
    group_by(region, lunaid) %>%
    filter(region %in% c('RCaudate','LCaudate','RThalamus')) %>%
    ungroup()
  
  luna_data <- luna_data %>%
    filter(!lunaid %in% c(11690, 11665))
  
  luna_data <- luna_data %>% group_by(lunaid,region) %>% mutate(nsess = 1:n()) %>% ungroup() %>% filter(nsess == 1)
  
  # wrangling
  df <- df %>% dplyr::select(RECID,hemi,doi_m,timepoint,age,GMrat,GPC,sex,roi,Glc,Glu,GPC.Cho,GABA,NAA,mI,Gln,NAAG,POSSX,Remitter_Status,Glu.Gln)
  luna_data <- luna_data %>% dplyr::select(lunaid,region,age,GMrat,GABA0,Glc0,Gln0,Glu0,mI0,NAA0,NAAG0,GPC.Cho0,GPC0,Glu.Gln0,sex)
  luna_data <- luna_data %>% mutate(group = 'HC', POSSX = NA, timepoint = 'BL')
  df <- df %>% dplyr::rename(id = 'RECID',GPC0 = 'GPC',Glu0 = 'Glu',GPC.Cho0 = 'GPC.Cho',Glc0 = 'Glc',GABA0 = 'GABA',NAA0 = 'NAA',mI0 = 'mI',Gln0 = 'Gln',NAAG0 = 'NAAG',Glu.Gln0 = 'Glu.Gln')
  luna_data <- luna_data %>% mutate(hemi = str_sub(region,1,1),
                                    roi = str_sub(region,2))
  luna_data <- luna_data %>% dplyr::select(!region)
  df <- df %>% dplyr::rename(group = "Remitter_Status")
  
  # merge
  luna_data <- luna_data %>% mutate(doi_m = NA) %>% rename(id = lunaid)
  df <- rbind(df,luna_data)
  df <- df %>% filter(roi %in% c('Caudate','Thalamus'))
  
  df$group <- relevel(factor(df$group), ref = "HC")
  df$timepoint <- relevel(factor(df$timepoint), ref = "BL")
  df$sex <- relevel(factor(df$sex),ref = "F")
  df <- df %>% mutate(age_sc = scale(age))
  
  df <- df %>%
    mutate(
      # Combine Group and Time into 3 valid conditions
      condition = case_when(
        group == "HC" ~ "HC_BL",
        group == "NR" & timepoint == "BL"  ~ "NR_BL",
        group == "NR" & timepoint == "FU"  ~ "NR_FU",
        group == "R" & timepoint == "BL"  ~ "R_BL",
        group == "R" & timepoint == "FU"  ~ "R_FU",      
      ),
      condition = factor(condition, levels = c("HC_BL", "NR_BL", "R_BL", "NR_FU","R_FU"))
    )
  
  df$condition <- relevel(df$condition, ref = "HC_BL")
  
  df <- df %>% mutate(group_level = case_when(group == 'HC' ~ 'HC',
                                              group == 'NR' ~ 'SZ',
                                              group == 'R' ~ 'SZ',
                                              is.na(group) ~ 'SZ')
  )
  df$group_level <- relevel(as.factor(df$group_level),'HC')
  
  
  
  lthal <- read_csv('L_Thalamus_HC_gamadj_long.csv')
  lthal <- lthal %>% rename(id = ld8) %>%
    mutate(doi_m = NA, hemi = 'L', timepoint = 'BL',roi = "Thalamus") %>% select("roi","doi_m","timepoint","id","hemi","age","GMrat","sex","Glu.Cr","Gln.Cr","GABA.Cr","GPC.Cr","NAA.Cr","NAAG.Cr","mI.Cr","Glu.Gln.Cr")
  lthal <- lthal %>% mutate(group = 'HC',condition = "HC_BL",group_level = "HC", POSSX = NA, GPC.Cho0 = NA, Glc0 = NA)
  lthal <- lthal %>% rename(GPC0 = "GPC.Cr",Glu0 = "Glu.Cr",Gln0 = "Gln.Cr",mI0 = "mI.Cr",GABA0 = "GABA.Cr",NAAG0 = "NAAG.Cr",NAA0 = "NAA.Cr",Glu.Gln0 = "Glu.Gln.Cr")
  lthal <- lthal %>% mutate(id = str_replace(id, "_[^_]*$", ""))
  
  df <- df %>% select(!age_sc)
  df <- rbind(df,lthal)  
  df <- df %>% mutate(age_sc = scale(age))
  df0 <- df %>% group_by(id) %>% slice(1) %>% ungroup()
  mDOI <- median(df0$doi_m[df0$group_level == 'SZ'],na.rm=TRUE)
  df <- df %>% mutate(median_split_doi = case_when(doi_m >= mDOI ~ 'high',
                                                   doi_m < mDOI ~ 'low')
  )
  df$median_split_doi[df$id == "2695"] = 'low'
  
  
  rthal <- read_csv('New_R_Thalamus_gamadj_long.csv') %>% select(!GABA.Cr_gamadj & !Gln.Cr_gamadj & !Glu.Cr_gamadj & !GPC.Cr_gamadj & !mI.Cr_gamadj & !NAA.Cr_gamadj & !NAAG.Cr_gamadj & !Glu.Gln.Cr_gamadj) %>% rename(GPC0 = "GPC.Cr",Glu0 = "Glu.Cr",Gln0 = "Gln.Cr",mI0 = "mI.Cr",GABA0 = "GABA.Cr",NAAG0 = "NAAG.Cr",NAA0 = "NAA.Cr",Glu.Gln0 = "Glu.Gln.Cr")
  rthal <- rthal %>% rename(id = RECID) %>%
    mutate(hemi = 'R', timepoint = 'BL',roi = "Thalamus") %>% select("age","GMrat","roi","timepoint","id","hemi","Glu.Gln0","NAA0","NAAG0","GABA0","mI0","Gln0","Glu0","GPC0")
  #rthal <- rthal %>% rename(SD_GABA = GABA.SD,SD_Glc = Glc.SD,SD_Gln = Gln.SD,SD_Glu = Glu.SD,SD_mI = mI.SD,SD_NAA = NAA.SD,SD_NAAG = NAAG.SD,SD_GPC.Cho = GPC.Cho.SD,SD_GPC = GPC.SD,SD_Glu.Gln = Glu.Gln.SD)
  #rthal <- rthal %>% select(!SD_Cre) #rename(SD_Cre = "Cre.SD")
  rthal <- rthal %>% mutate(id = str_replace(id, "_[^_]*$", ""))
  rthal <- rthal %>% mutate(doi_m = NA, median_split_doi = NA, sex = NA, GPC.Cho0 = NA, Glc0 = NA, POSSX = NA, group = 'HC', condition = 'BL',group_level = 'HC', age_sc = NA)
  
  df <- df %>% filter(!(roi == 'Thalamus' & hemi == 'R' & group_level == 'HC'))
  df <- df %>% group_by(id,roi,hemi,timepoint) %>% slice(1) %>% ungroup()
  rthal <- rthal %>% group_by(id,roi,hemi,timepoint) %>% slice(1) %>% ungroup()
  df <- rbind(df,rthal)
  
  df <- df %>% group_by(group_level,id,timepoint,roi,hemi) %>% slice(1) %>% ungroup()
  
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu0 ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA0 ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), mI0 ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Gln0 ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln0 ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG0 ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAA0 ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus'), GPC0 ~ group_level*hemi + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group_level + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGln,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  
  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu0~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA0 ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), mI0 ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Gln0 ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu.Gln0 ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAAG0 ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAA0 ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate'), GPC0 ~ group_level*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GPC.Cho ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGln,CaGluGln,CaNAAG,CaNAA,CaGpc) %>% mutate(roi = 'Caudate')
  
  Mall_R_SZ <- rbind(Mth1,Mca1) %>% filter(term %in% c('group_levelSZ:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr)
  Mall_SZ <- rbind(Mth1,Mca1) %>% filter(term %in% c('group_levelSZ')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr)
  
  
}

if (reload==T){
  
  #read in ALL corrected data (including duplicates)
  #isolate patients with baseline and follow-up only (ie only patients with known remitter status)
  sarpal_data_all <- read_csv('sarpal_mrsi_adj_all_20250310.csv') %>%
    filter(source %in% c('BL.FU','DARES.baselines')) %>% group_by(RECID,timepoint,hemi,roi) %>% slice(1) %>% ungroup()
  
  #add in sex from outside spreadsheet
  # 2026-06-02 AndyP changed to sheet 3 was reading incorrect sheet
  sarpal_data_supplemental <- read_excel('7T.DARES.baselines.xlsx', sheet = 3) %>% group_by(RECID) %>% slice(1) %>% ungroup()
  
  df <- inner_join(sarpal_data_all,sarpal_data_supplemental,by='RECID')
  
  df <- df %>% mutate(sex = as.factor(case_when(sex == 1 ~ 'M',
                                                sex == 2 ~ 'F')))
  
  
  ############################
  #add in remitter status from outside spreadsheet
  supplemental_data2 <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_data2 <- supplemental_data2 %>% filter(Timepoint != 0)
  supplemental_data2 <- supplemental_data2 %>%
    mutate(Remitter_Status = ifelse(CONCDIS4 <= 3 & HALL12 <= 3 & UTC15 <= 3, "R", "NR"), timepoint = Timepoint)
  supplemental_data2 <- supplemental_data2 %>% select(RECID,POSSX,Remitter_Status, timepoint)
  supplemental_dataBL <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_dataBL <- supplemental_dataBL %>% filter(Timepoint == 0) %>% mutate(timepoint = Timepoint, Remitter_Status = NA)
  supplemental_dataBL <- supplemental_dataBL %>% select(RECID, POSSX,timepoint,Remitter_Status)
  supplemental_data2 <- rbind(supplemental_data2, supplemental_dataBL)
  supplemental_data2 <- supplemental_data2 %>% group_by(RECID) %>% tidyr::fill(Remitter_Status, .direction = "downup") %>% ungroup() %>% mutate(timepoint = case_when(timepoint == 0 ~ 1,
                                                                                                                                                                      timepoint > 0 ~ 2))
  # sarpal_data_BLFU <- sarpal_data_BLFU %>% 
  #   mutate(Remitter_Status = supplemental_data2$Remitter_Status[match(RECID,supplemental_data2$RECID)])
  
  df <- left_join(df,supplemental_data2,by=c('RECID','timepoint'))
  df <- df %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                            timepoint == 2 ~ 'FU'))
  
  df$doi_m <- df$`DUP (months)`
  # note Cr_gamadj is Metabolite / Creatine NOT raw Creatine metabolite.
  # The adjustment is based on the length of the scan and was done for luna lab and Sarpal lab.
  load('luna_data.Rdata')
  #get rid of unimportant columns for our immediate analysis to allow for wide pivot
  luna_data <- luna_data[, c("met", "ld8", "region", "age", "GMrat","Cr_gamadj")]
  #pivot wide for easy comparison with other spreadsheet
  luna_data <- luna_data %>%
    pivot_wider(names_from = met, values_from = Cr_gamadj)
  luna_data <- luna_data %>%
    separate(ld8, into = c("lunaid", "date"), sep = "_")
  luna_data$lunaid <- as.numeric(luna_data$lunaid)
  luna_data <- luna_data %>% filter(age >= 18)
  
  #add in sex from outside spreadsheet
  luna_data_supplemental <- read_csv('luna7t_sifpc_subset.csv')
  luna_data_supplemental <- luna_data_supplemental %>% dplyr::select(lunaid,sex) %>% group_by(lunaid) %>% slice(1)
  
  luna_data <- inner_join(luna_data,luna_data_supplemental,by=c("lunaid"))
  
  luna_data <- luna_data %>%
    mutate(date = as.numeric(date)) %>%
    arrange(region, lunaid, date) %>%
    group_by(region, lunaid) %>%
    filter(region %in% c('RCaudate','LCaudate','RThalamus')) %>%
    ungroup()
  
  luna_data <- luna_data %>%
    filter(!lunaid %in% c(11690, 11665))
  
  luna_data <- luna_data %>% group_by(lunaid,region) %>% mutate(nsess = 1:n()) %>% ungroup() %>% filter(nsess == 1)
  
  # wrangling
  df <- df %>% dplyr::select(RECID,hemi,doi_m,timepoint,age,GMrat,GPC.adj,sex,roi,Glc.adj,Glu.adj,GPC.Cho.adj,GABA.adj,NAA.adj,mI.adj,Gln.adj,NAAG.adj,POSSX,Remitter_Status,Glu.Gln.adj)
  luna_data <- luna_data %>% dplyr::select(lunaid,date,region,age,GMrat,GABA,Glc,Gln,Glu,mI,NAA,NAAG,GPC.Cho,GPC,Glu.Gln,sex)
  luna_data <- luna_data %>% mutate(group = 'HC', POSSX = NA, timepoint = 'BL')
  df <- df %>% dplyr::rename(id = 'RECID',GPC = 'GPC.adj',Glu = 'Glu.adj',GPC.Cho = 'GPC.Cho.adj',Glc = 'Glc.adj',GABA = 'GABA.adj',NAA = 'NAA.adj',mI = 'mI.adj',Gln = 'Gln.adj',NAAG = 'NAAG.adj',Glu.Gln = 'Glu.Gln.adj')
  luna_data <- luna_data %>% mutate(hemi = str_sub(region,1,1),
                                    roi = str_sub(region,2))
  luna_data <- luna_data %>% dplyr::select(!region)
  df <- df %>% dplyr::rename(group = "Remitter_Status")
  # merge
  luna_data <- luna_data %>% mutate(doi_m = NA) %>% rename(id = lunaid)
  df$date <- NA
  df <- rbind(df,luna_data)
  df <- df %>% filter(roi %in% c('Caudate','Thalamus'))
  
  df$group <- relevel(factor(df$group), ref = "HC")
  df$timepoint <- relevel(factor(df$timepoint), ref = "BL")
  df$sex <- relevel(factor(df$sex),ref = "F")
  df <- df %>% mutate(age_sc = scale(age))
  
  df <- df %>%
    mutate(
      # Combine Group and Time into 3 valid conditions
      condition = case_when(
        group == "HC" ~ "HC_BL",
        group == "NR" & timepoint == "BL"  ~ "NR_BL",
        group == "NR" & timepoint == "FU"  ~ "NR_FU",
        group == "R" & timepoint == "BL"  ~ "R_BL",
        group == "R" & timepoint == "FU"  ~ "R_FU",      
      ),
      condition = factor(condition, levels = c("HC_BL", "NR_BL", "R_BL", "NR_FU","R_FU"))
    )
  
  df$condition <- relevel(df$condition, ref = "HC_BL")
  
  df <- df %>% mutate(group_level = case_when(group == 'HC' ~ 'HC',
                                              group == 'NR' ~ 'SZ',
                                              group == 'R' ~ 'SZ',
                                              is.na(group) ~ 'SZ')
                      )
  df$group_level <- relevel(as.factor(df$group_level),'HC')
  
  table <- df %>% group_by(roi,hemi,group_level,timepoint) %>% summarize(mgaba = mean(GABA,na.rm=TRUE), 
                                                         sgaba = sd(GABA,na.rm=TRUE),
                                                         mglu = mean(Glu, na.rm=TRUE),
                                                         sglu = sd(Glu, na.rm=TRUE),
                                                         mnaa = mean(NAA,na.rm=TRUE),
                                                         snaa = sd(NAA,na.rm=TRUE),
                                                         mmI = mean(mI,na.rm=TRUE),
                                                         smI = sd(mI,na.rm=TRUE),
                                                         mnaag = mean(NAAG,na.rm=TRUE),
                                                         snaag = sd(NAAG,na.rm=TRUE)
  ) %>% ungroup() %>% arrange(group_level,hemi)
  
  
  lthal <- read_csv('L_Thalamus_HC_gamadj_long.csv')
  lthal <- lthal %>% rename(id = ld8) %>%
    mutate(doi_m = NA, hemi = 'L', timepoint = 'BL',roi = "Thalamus") %>% select("roi","doi_m","timepoint","id","hemi","age","GMrat","sex","Glu.Cr_gamadj","Gln.Cr_gamadj","GABA.Cr_gamadj","GPC.Cr_gamadj","NAA.Cr_gamadj","NAAG.Cr_gamadj","mI.Cr_gamadj","Glu.Gln.Cr_gamadj")
  lthal <- lthal %>% mutate(group = 'HC',condition = "HC_BL",group_level = "HC", POSSX = NA, GPC.Cho = NA, Glc = NA)
  lthal <- lthal %>% rename(GPC = "GPC.Cr_gamadj",Glu = "Glu.Cr_gamadj",Gln = "Gln.Cr_gamadj",mI = "mI.Cr_gamadj",GABA = "GABA.Cr_gamadj",NAAG = "NAAG.Cr_gamadj",NAA = "NAA.Cr_gamadj",Glu.Gln = "Glu.Gln.Cr_gamadj")
  lthal <- lthal %>% mutate(id = str_replace(id, "_[^_]*$", ""))
  lthal$date <- NA
  df <- df %>% select(!age_sc)
  df <- rbind(df,lthal)  
  df <- df %>% mutate(age_sc = scale(age))
  df0 <- df %>% group_by(id) %>% slice(1) %>% ungroup()
  mDOI <- median(df0$doi_m[df0$group_level == 'SZ'],na.rm=TRUE)
  df <- df %>% mutate(median_split_doi = case_when(doi_m >= mDOI ~ 'high',
                                                   doi_m < mDOI ~ 'low')
  )
  df$median_split_doi[df$id == "2695"] = 'low'
  #cre1 <- read_csv('Cr_raw_all_gamadj_long.csv') # missing left thalamus
  cre <- read_csv('metabolites_gamadj_SSD_long.csv')
  cre <- cre %>% rename(id = RECID) %>%
    mutate(hemi = case_when(region == 'left caudate' ~ 'L',
                            region == 'right caudate' ~ 'R',
                            region == 'right thalamus' ~ 'R',
                            region == 'left thalamus' ~ 'L')) %>%
    mutate(roi = case_when(biregion == 'thalamus' ~ 'Thalamus',
                           biregion == 'caudate' ~ 'Caudate')) %>%
    mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                 timepoint == 2 ~ 'FU')) %>%
    select("roi","timepoint","id","hemi",'Cre_gamadj',"Cre..SD") %>%
    rename(Cre = "Cre_gamadj", SD_Cre = "Cre..SD") 
  
  cre_hc <- read_csv('Cr_raw_gamadj_long_HC.csv')
  cre_hc <- cre_hc %>% rename(id = ld8) %>%
    mutate(hemi = case_when(region == 'LCaudate' ~ 'L',
                            region == 'RCaudate' ~ 'R',
                            region == 'RThalamus' ~ 'R')) %>%
    mutate(roi = case_when(biregion == 'Thalamus' ~ 'Thalamus',
                           biregion == 'Caudate' ~ 'Caudate')) %>%
    group_by(id,roi,hemi) %>% mutate(SES = n()) %>% filter(SES == 1) %>%
    mutate(timepoint = case_when(SES == '1' ~ 'BL')) %>%
    select("roi","timepoint","id","hemi",'Cr_raw_gamadj','Cr_raw.SD') %>%
    rename(Cre = "Cr_raw_gamadj", SD_Cre = "Cr_raw.SD") %>%
    filter(!(roi == 'Thalamus' & hemi == 'R'))
  
  # reloading to get Cre and CRLB (SD)
  lthal <- read_csv('L_Thalamus_HC_gamadj_long.csv')
  lthal <- lthal %>% rename(id = ld8) %>%
    mutate(hemi = 'L', timepoint = 'BL',roi = "Thalamus") %>% select("roi","timepoint","id","hemi","Cr_raw_gamadj","Cre.SD","GABA.SD","Glc.SD","Gln.SD","Glu.SD","mI.SD","NAA.SD","NAAG.SD","GPC.SD","Glu.Gln.SD","GPC.Cho.SD")
  lthal <- lthal %>% rename(SD_GABA = GABA.SD,SD_Glc = Glc.SD,SD_Gln = Gln.SD,SD_Glu = Glu.SD,SD_mI = mI.SD,SD_NAA = NAA.SD,SD_NAAG = NAAG.SD,SD_GPC.Cho = GPC.Cho.SD,SD_GPC = GPC.SD,SD_Glu.Gln = Glu.Gln.SD)
  lthal <- lthal %>% rename(Cre = "Cr_raw_gamadj",SD_Cre = "Cre.SD")
  lthal <- lthal %>% mutate(id = str_replace(id, "_[^_]*$", ""))
  
  rthal <- read_csv('New_R_Thalamus_gamadj_long.csv') %>% select(!GABA & !Gln & !Glu & !GPC & !mI & !NAA & !NAAG & !Glu.Gln & !Cre) %>% rename(GPC = "GPC.Cr_gamadj",Glu = "Glu.Cr_gamadj",Gln = "Gln.Cr_gamadj",mI = "mI.Cr_gamadj",GABA = "GABA.Cr_gamadj",NAAG = "NAAG.Cr_gamadj",NAA = "NAA.Cr_gamadj",Glu.Gln = "Glu.Gln.Cr_gamadj")
  rthal <- rthal %>% rename(id = RECID) %>%
    mutate(hemi = 'R', timepoint = 'BL',roi = "Thalamus") %>% select("age","GMrat","roi","timepoint","id","hemi","Glu.Gln","NAA","NAAG","GABA","mI","Gln","Glu","GPC")
  #rthal <- rthal %>% rename(SD_GABA = GABA.SD,SD_Glc = Glc.SD,SD_Gln = Gln.SD,SD_Glu = Glu.SD,SD_mI = mI.SD,SD_NAA = NAA.SD,SD_NAAG = NAAG.SD,SD_GPC.Cho = GPC.Cho.SD,SD_GPC = GPC.SD,SD_Glu.Gln = Glu.Gln.SD)
  #rthal <- rthal %>% select(!SD_Cre) #rename(SD_Cre = "Cre.SD")
  rthal <- rthal %>% mutate(id = str_replace(id, "_[^_]*$", ""))
  rthal <- rthal %>% mutate(doi_m = NA, median_split_doi = NA, sex = NA, GPC.Cho = NA, Glc = NA, POSSX = NA, group = 'HC', condition = 'BL',group_level = 'HC', age_sc = NA)
  
  df <- df %>% filter(!(roi == 'Thalamus' & hemi == 'R' & group_level == 'HC'))
  df <- df %>% group_by(id,roi,hemi,timepoint) %>% slice(1) %>% ungroup()
  rthal <- rthal %>% group_by(id,roi,hemi,timepoint) %>% slice(1) %>% ungroup()
  rthal$date <- NA
  df <- rbind(df,rthal)
  
  
  rthal <- read_csv('New_R_Thalamus_gamadj_long.csv')
  rthal <- rthal %>% rename(id = RECID) %>%
    mutate(hemi = 'R', timepoint = 'BL',roi = "Thalamus") %>% select("roi","timepoint","id","hemi","Cre","Cre.SD","GABA.SD","Glc.SD","Gln.SD","Glu.SD","mI.SD","NAA.SD","NAAG.SD","GPC.SD","Glu.Gln.SD","GPC.Cho.SD")
  rthal <- rthal %>% rename(SD_GABA = GABA.SD,SD_Glc = Glc.SD,SD_Gln = Gln.SD,SD_Glu = Glu.SD,SD_mI = mI.SD,SD_NAA = NAA.SD,SD_NAAG = NAAG.SD,SD_GPC.Cho = GPC.Cho.SD,SD_GPC = GPC.SD,SD_Glu.Gln = Glu.Gln.SD)
  rthal <- rthal %>% rename(SD_Cre = "Cre.SD")
  rthal <- rthal %>% mutate(id = str_replace(id, "_[^_]*$", ""))
  
  cre <- rbind(cre,cre_hc)
  cre <- cre %>% mutate(SD_GABA = NA,SD_Glc = NA,SD_Gln = NA,SD_Glu = NA,SD_mI = NA,SD_NAA = NA,SD_NAAG = NA,SD_GPC.Cho = NA,SD_GPC = NA,SD_Glu.Gln = NA)
  cre <- rbind(cre,lthal)
  cre <- rbind(cre,rthal)
  cre <- cre %>% mutate(group_level = case_when(as.numeric(id) < 10000 ~ 'SZ',
                                                as.numeric(id) >= 10000 ~ 'HC')
  )
  #hc_exclude <- c("11390_20180628","11632_20180416","11665_20180628","11701_20181015","11702_20181109","11754_20210823","11788_20200117")
  #df <- df %>% filter(!id %in% hc_exclude)
  cre$id <- as.character(cre$id)
  cre <- cre %>% arrange(id,timepoint,roi,hemi)
  
  
  # handedness
  hand <- data.frame(id = c("2207","2229","2245","2246","2276","2279","2292","2318","2326","2343","2346","2351","2358","2367","2379","2497","2578","2647","2690","2773","2824","2859"), hand = c('L','R','R','R','L','R','R','R','M','R','R','L','R',NA,NA,'R','R','R',NA,'R',NA,NA))
  
  # testing L vs R CRLB and metabolite ratios (just for testing).  For some reason, Caudate is not getting picked up correctly in this merge.
  test <- full_join(df,cre,by=c('id','group_level','timepoint','hemi','roi'))
  test <- test %>% group_by(id,hemi,roi,timepoint) %>% slice(1) %>% ungroup()
  df <- left_join(df,hand,by=c('id'))
  
  Th <- test %>% filter(roi == 'Thalamus') %>% ungroup()
  Th <- Th %>% select(id,hemi,roi,NAAG,GABA,NAA,Cre,GMrat,group_level,timepoint,SD_NAAG,SD_GABA,SD_NAA,SD_Cre) %>% 
    pivot_wider(id_cols = c('id','group_level'),names_from = c('timepoint','hemi'), values_from = c('NAAG','NAA','GABA','Cre','GMrat','SD_NAAG','SD_GABA','SD_NAA','SD_Cre'))
  ggplot(data=Th) + 
    geom_point(aes(x=NAAG_BL_R, y=NAAG_BL_L,color = group_level)) + 
    geom_point(aes(x=NAAG_FU_R, y=NAAG_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right NAAG') + ylab('Left NAAG') + 
    geom_abline(intercept = 0, slope = 1) + xlim(c(0,1.3)) + ylim(c(0,1.5))
  ggplot(data=Th) + 
    geom_point(aes(x=NAA_BL_R, y=NAA_BL_L,color = group_level)) + 
    geom_point(aes(x=NAA_FU_R, y=NAA_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right NAA') + ylab('Left NAA') + 
    geom_abline(intercept = 0, slope = 1) + xlim(c(0.75,1.75)) + ylim(c(0.75,1.75)) 
  ggplot(data=Th) + 
    geom_point(aes(x=Cre_BL_R, y=Cre_BL_L,color = group_level)) + 
    geom_point(aes(x=Cre_FU_R, y=Cre_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right Cre') + ylab('Left Cre') + 
    geom_abline(intercept = 0, slope = 1)   
  
  
  ggplot(data=Th,aes(label = id)) + 
    geom_point(aes(x=GMrat_BL_R, y=GMrat_BL_L,color = group_level)) + 
    geom_point(aes(x=GMrat_FU_R, y=GMrat_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right GMrat') + ylab('Left GMrat') + 
    geom_abline(intercept = 0, slope = 1) + geom_text_repel(data=Th,aes(x=GMrat_BL_R, y=GMrat_BL_L),size = 3, box.padding = 0.1) 
  
  
  # plot CRLB
  
  ggplot(data=Th, aes(label = id)) + 
    geom_jitter(aes(x=SD_Cre_BL_R, y=SD_Cre_BL_L,color = group_level)) + 
    geom_jitter(aes(x=SD_Cre_FU_R, y=SD_Cre_FU_R, color = 'follow-up')) + 
    theme(aspect.ratio = 1) + geom_text_repel(data=Th,aes(x=SD_Cre_BL_R, y=SD_Cre_BL_L),size = 3, box.padding = 0.1) +
    xlab('Right SD Cre Thalamus') + ylab('Left SD Cre Thalamus') + 
    geom_abline(intercept = 0, slope = 1) + xlim(c(0,25)) + ylim(c(0,25))
  
  ggplot(data=Th, aes(label = id)) + 
    geom_jitter(aes(x=SD_NAAG_BL_R, y=SD_NAAG_BL_L,color = group_level)) + 
    geom_jitter(aes(x=SD_NAAG_FU_R, y=SD_NAAG_FU_R, color = 'follow-up')) + 
    geom_text_repel(data=Th,aes(x=SD_NAAG_BL_R, y=SD_NAAG_BL_L),size = 3, box.padding = 0.1) +   
    theme(aspect.ratio = 1) + 
    xlab('Right SD NAAG Thalamus') + ylab('Left SD NAAG Thalamus') + 
    geom_abline(intercept = 0, slope = 1)  + xlim(c(0,25)) + ylim(c(0,25))
  
  ggplot(data=Th, aes(label = id)) + 
    geom_jitter(aes(x=SD_NAA_BL_R, y=SD_NAA_BL_L,color = group_level)) + 
    geom_jitter(aes(x=SD_NAA_FU_R, y=SD_NAA_FU_R, color = 'follow-up')) + 
    geom_text_repel(data=Th,aes(x=SD_NAA_BL_R, y=SD_NAA_BL_L),size = 3, box.padding = 0.1) +
    theme(aspect.ratio = 1) + 
    xlab('Right SD NAA Thalamus') + ylab('Left SD NAA Thalamus') + 
    geom_abline(intercept = 0, slope = 1) + xlim(c(0,25)) + ylim(c(0,25)) 
  
  Ca <- test %>% filter(roi == 'Caudate')
  Ca <- Ca %>% select(id,hemi,roi,NAAG,GABA,NAA,Cre,GMrat,group_level,timepoint,SD_GABA,SD_NAA,SD_Cre) %>%
    pivot_wider(id_cols = c('id','group_level'),names_from = c('timepoint','hemi'), values_from = c('NAAG','NAA','GABA','Cre','GMrat','SD_GABA','SD_NAA','SD_Cre'))
  ggplot(data=Ca) +
    geom_point(aes(x=NAAG_BL_R, y=NAAG_BL_L,color = group_level)) +
    geom_point(aes(x=NAAG_FU_R, y=NAAG_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right NAAG') + ylab('Left NAAG') +
    geom_abline(intercept = 0, slope = 1) + xlim(c(0,1.3)) + ylim(c(0,1.5))
  ggplot(data=Ca) +
    geom_point(aes(x=NAA_BL_R, y=NAA_BL_L,color = group_level)) +
    geom_point(aes(x=NAA_FU_R, y=NAA_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right NAA') + ylab('Left NAA') +
    geom_abline(intercept = 0, slope = 1) + xlim(c(0.75,1.75)) + ylim(c(0.75,1.75))
  ggplot(data=Ca) +
    geom_point(aes(x=Cre_BL_R, y=Cre_BL_L,color = group_level)) +
    geom_point(aes(x=Cre_FU_R, y=Cre_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right Cre') + ylab('Left Cre') +
    geom_abline(intercept = 0, slope = 1)
  ggplot(data=Ca) +
    geom_point(aes(x=GMrat_BL_R, y=GMrat_BL_L,color = group_level)) +
    geom_point(aes(x=GMrat_FU_R, y=GMrat_FU_L, color = 'follow-up')) + theme(aspect.ratio = 1) + xlab('Right GMrat') + ylab('Left GMrat') +
    geom_abline(intercept = 0, slope = 1)

  ggplot(data=Ca, aes(label = id)) +
    geom_jitter(aes(x=SD_Cre_BL_R, y=SD_Cre_BL_L,color = group_level)) +
    geom_jitter(aes(x=SD_Cre_FU_R, y=SD_Cre_FU_R, color = 'follow-up')) +
    geom_text_repel(data=Ca,aes(x=SD_Cre_BL_R, y=SD_Cre_BL_L),size = 3, box.padding = 0.1) +
    theme(aspect.ratio = 1) +
    xlab('Right SD Cre Caudate') + ylab('Left SD Cre Caudate') +
    geom_abline(intercept = 0, slope = 1)  + xlim(c(0,25)) + ylim(c(0,25))

  ggplot(data=Ca,aes(label = id)) +
    geom_jitter(aes(x=SD_NAA_BL_R, y=SD_NAA_BL_L,color = group_level)) +
    geom_jitter(aes(x=SD_NAA_FU_R, y=SD_NAA_FU_R, color = 'follow-up')) +
    geom_text_repel(data=Ca,aes(x=SD_NAA_BL_R, y=SD_NAA_BL_L),size = 3, box.padding = 0.1) +
    theme(aspect.ratio = 1) +
    xlab('Right SD NAA Caudate') + ylab('Left SD NAA Caudate') +
    geom_abline(intercept = 0, slope = 1)
  
}
#####################################
### Group Level R/NR by time.     ###
#####################################

if (longitudinal==T){
  
  df <- df %>% filter(group_level == 'SZ' | (group_level == 'HC' & timepoint == 'BL'))
  
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGlu0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), Glu ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThGABA0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), GABA ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), mI ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThmI0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), mI ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  
  #ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  #ThGln0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), Gln ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThGluGln0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), Glu.Gln ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAA ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  ThNAA0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), NAA ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect)  %>% mutate(metabolite = 'NAA')
  
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAAG0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), NAAG ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect)  %>% mutate(metabolite = 'NAAG')
  
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus'), GPC ~ condition*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  ThGpc0 <- tidy(lm(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), GPC ~ group*timepoint + sex + scale(GMrat)))  %>% mutate(metabolite = 'GPC')
  
  # ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ condition + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  # ThCho0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), GPC.Cho ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect)  %>% mutate(metabolite = 'GPC.Cho')
  # 
  Mth <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  #############################
  ####### Caudate #############
  #############################
  
  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGlu0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), Glu ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CaGABA0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), GABA ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), mI ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CamI0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), mI ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  
  #CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  #CaGln0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), Gln ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  
  CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu.Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaGluGln0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), Glu.Gln ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  
  CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAAG ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAAG0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), NAAG ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect)  %>% mutate(metabolite = 'NAAG')
  
  CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAA ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  CaNAA0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), NAA ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect)  %>% mutate(metabolite = 'NAA')
  
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate'), GPC ~ condition*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaGpc0 <- tidy(lm(data = df %>% filter(group != 'HC' & roi == 'Caudate'), GPC ~ group*timepoint + sex + scale(GMrat)))  %>% mutate(metabolite = 'GPC')
  
  CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GPC.Cho ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  CaCho0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), GPC.Cho ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect)  %>% mutate(metabolite = 'GPC.Cho')
  
  Mca <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAAG,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  
  # 9 metabolites x 2 regions = N/18 for FDR
  Mall_NR_BL <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_BL')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_R_BL <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_BL')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_NR_FU <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_FU')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_R_FU <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_FU')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_NR_BL_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_BL:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_R_BL_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_BL:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_NR_FU_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_FU:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_R_FU_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_FU:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  
  
  M <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAA ~ condition*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ condition | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  M <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ condition | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  M <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu.Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ condition | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  M <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ condition*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ condition | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
}
############################
#### High/Low DOI ##########
############################

if (hilowdoi==T){

  df <- df %>% mutate(doigr = case_when(group_level == 'HC' ~ 'HC',
                                        group_level == 'SZ' ~ median_split_doi))
  df$doigr <- relevel(as.factor(df$doigr),ref='HC')                      

  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAA ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  #ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GPC ~ doigr*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')


  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  #CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAAG ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate'), GPC ~ doigr*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC.Cho ~ doigr*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')

  Mall_high_r <- rbind(Mth1,Mca1) %>% filter(term %in% c('doigrhigh:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_high <- rbind(Mth1,Mca1) %>% filter(term %in% c('doigrhigh')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_low_r <- rbind(Mth1,Mca1) %>% filter(term %in% c('doigrlow:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_low <- rbind(Mth1,Mca1) %>% filter(term %in% c('doigrlow')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)

  M <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'),  Glu ~ doigr*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ doigr | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  
  Th <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ doigr*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(Th, ~ doigr | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
}


#############################
### Group Level ###
#############################

if (group==T){
  
  
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Gln ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAA ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  #ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GPC ~ group_level*hemi + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group_level + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')
  
  
  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  #CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC ~ group_level*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC.Cho ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  Mall_R_SZ <- rbind(Mth1,Mca1) %>% filter(term %in% c('group_levelSZ:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr) %>% filter(pfdr < 0.05)
  Mall_SZ <- rbind(Mth1,Mca1) %>% filter(term %in% c('group_levelSZ')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr) %>% filter(pfdr < 0.05)
  
  
  library(emmeans)
  
  # examine emmeans
  MCaGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'),  Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGlu, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")

  # fv <- fitted(MCaGlu)
  # rsq <- (residuals(MCaGlu))^2
  # valid_indices <- which(fv > 0 & rsq > 0)
  # ln_res_sq <- log(rsq[valid_indices])
  # ln_fitted <- log(fv[valid_indices])
  # park_model <- lm(ln_res_sq ~ ln_fitted)
  # summary(park_model)
  # 
  # MCaGlu <- glmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'),  Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id), family = inverse.gaussian(link = "1/mu^2"))
  # emm <- emmeans(MCaGlu, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # # Convert to data frame
  # emm_df <- as.data.frame(emm)
  # pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCaGPC.Cho <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'),  GPC.Cho ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGPC.Cho, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
    
  # examine emmeans
  MThGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  fv <- fitted(MThGABA)
  rsq <- (residuals(MThGABA))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  # examine emmeans
  MThGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGlu, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  fv <- fitted(MThGlu)
  rsq <- (residuals(MThGlu))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGluGln, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
 
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGluGln, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGluGln, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThmI <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThmI, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # this one is -2 => Variance is inversely proportional to the mean squared, inverse heteroscedasticity
  fv <- fitted(MThmI)
  rsq <- (residuals(MThmI))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  # How to fix:
  # 1. Extract the fitted gamma from your Park model
  # gamma_hat <- coef(park_model)["ln_fitted"]
  # 
  # # 2. Calculate the weights (Inverse of variance)
  # # Since Var = Mean^Gamma, Weight = 1 / (Mean^Gamma)
  # weights_wls <- 1 / (fitted_vals^gamma_hat)
  # 
  # # 3. Re-run your baseline model with weights
  # wls_model <- lm(y ~ x1 + x2, data = your_data, weights = weights_wls)
  # summary(wls_model)
  
  # examine emmeans
  MCamI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCamI, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # this one is significant, variance scales linearly with the mean, consider using Poisson GLM
  fv <- fitted(MCamI)
  rsq <- (residuals(MCamI))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # N.S.
  fv <- fitted(MThNAAG)
  rsq <- (residuals(MThNAAG))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  # examine emmeans, yes there is more in L in HC and more in R in SZ
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAAG, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAA, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # N.S.
  fv <- fitted(MThNAA)
  rsq <- (residuals(MThNAA))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
   
  # examine emmeans
  MCaNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaNAA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCaNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaNAA, ~ hemi | group_level, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCamI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCamI, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  # # Create a customized plot
  # ggplot(emm_df, aes(x = group_level, y = emmean, color = hemi, group=hemi)) +
  #   geom_point(size = 3) +
  #   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  #   theme_minimal() +
  #   labs(title = "Model-Predicted Means by Group",
  #        y = "Estimated Marginal Mean",
  #        x = "Group")
  MCaGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGluGln, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  fv <- fitted(MCaGluGln)
  rsq <- (residuals(MCaGluGln))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  MCaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
}

#####################
###.  E/I Balance ###
#####################

if (eibalance){
  
  # df <- df %>% group_by(group_level) %>% mutate(GABA_sc = scale(GABA), Glu_sc = scale(Glu))
  # 
  # 
  # #EI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ scale(Glu.Gln)*group_level*hemi + sex + scale(GMrat) + (1|id))
  # EI1 <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & group != 'HC'), GABA_sc ~ Glu_sc:group:timepoint + group*timepoint + hemi + age_sc + sex + scale(GMrat) + (1|id))
  # #EI2 <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ scale(Gln)*group_level*hemi + sex + scale(GMrat) + (1|id))
  # 
  # 
  # #EI <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ scale(Glu.Gln)*group + hemi + sex + scale(GMrat) + (1|id))
  # EI2 <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & group != 'HC'), GABA_sc ~ Glu_sc*group + timepoint + hemi + sex + scale(GMrat) + (1|id))
  # #EI2 <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ scale(Gln)*group + hemi + sex + scale(GMrat) + (1|id))
  # 
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_HC_Caudate_BL))
  b <- as.numeric(unlist(dfw$GABA_HC_Caudate_BL))
  c <- as.numeric(unlist(dfw$Glu_SZ_Caudate_BL))
  d <- as.numeric(unlist(dfw$GABA_SZ_Caudate_BL))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Caudate = a1, GABA_HC_Caudate = b1), 
              SZ = data.frame(Glu_SZ_Caudate = c1, GABA_SZ_Caudate = d1))
  
  ccCa <- cocor(~Glu_HC_Caudate + GABA_HC_Caudate | Glu_SZ_Caudate + GABA_SZ_Caudate, dfw, alternative = 'two.sided')
  print(ccCa)
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_HC_Thalamus_BL))
  b <- as.numeric(unlist(dfw$GABA_HC_Thalamus_BL))
  c <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_BL))
  d <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_BL))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Thalamus = a1, GABA_HC_Thalamus = b1), 
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))
  
  ccTh <- cocor(~Glu_HC_Thalamus + GABA_HC_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')
  print(ccTh)
  #####################
  #### left side #####
  #####################
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_HC_Caudate_BL_L))
  b <- as.numeric(unlist(dfw$GABA_HC_Caudate_BL_L))
  c <- as.numeric(unlist(dfw$Glu_SZ_Caudate_BL_L))
  d <- as.numeric(unlist(dfw$GABA_SZ_Caudate_BL_L))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Caudate_L = a1, GABA_HC_Caudate_L = b1), 
              SZ = data.frame(Glu_SZ_Caudate_L = c1, GABA_SZ_Caudate_L = d1))
  
  ccCa <- cocor(~Glu_HC_Caudate_L + GABA_HC_Caudate_L | Glu_SZ_Caudate_L + GABA_SZ_Caudate_L, dfw, alternative = 'two.sided')
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_HC_Thalamus_BL_L))
  b <- as.numeric(unlist(dfw$GABA_HC_Thalamus_BL_L))
  c <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_BL_L))
  d <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_BL_L))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Thalamus_L = a1, GABA_HC_Thalamus_L = b1), 
              SZ = data.frame(Glu_SZ_Thalamus_L = c1, GABA_SZ_Thalamus_L = d1))
  
  ccTh <- cocor(~Glu_HC_Thalamus_L + GABA_HC_Thalamus_L | Glu_SZ_Thalamus_L + GABA_SZ_Thalamus_L, dfw, alternative = 'two.sided')
  
  #####################
  #### right side #####
  #####################
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_HC_Caudate_BL_R))
  b <- as.numeric(unlist(dfw$GABA_HC_Caudate_BL_R))
  c <- as.numeric(unlist(dfw$Glu_SZ_Caudate_BL_R))
  d <- as.numeric(unlist(dfw$GABA_SZ_Caudate_BL_R))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Caudate_R = a1, GABA_HC_Caudate_R = b1), 
              SZ = data.frame(Glu_SZ_Caudate_R = c1, GABA_SZ_Caudate_R = d1))
  
  ccCa <- cocor(~Glu_HC_Caudate_R + GABA_HC_Caudate_R | Glu_SZ_Caudate_R + GABA_SZ_Caudate_R, dfw, alternative = 'two.sided')
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_HC_Thalamus_BL_R))
  b <- as.numeric(unlist(dfw$GABA_HC_Thalamus_BL_R))
  c <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_BL_R))
  d <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_BL_R))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Thalamus_R = a1, GABA_HC_Thalamus_R = b1), 
              SZ = data.frame(Glu_SZ_Thalamus_R = c1, GABA_SZ_Thalamus_R = d1))
  
  ccTh <- cocor(~Glu_HC_Thalamus_R + GABA_HC_Thalamus_R | Glu_SZ_Thalamus_R + GABA_SZ_Thalamus_R, dfw, alternative = 'two.sided')
  
  
  # ########################
  # ### timepoint ##########
  # ########################
  # 
  # 
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)

  a <- as.numeric(unlist(dfw$Glu_SZ_Caudate_BL))
  b <- as.numeric(unlist(dfw$GABA_SZ_Caudate_BL))
  c <- as.numeric(unlist(dfw$Glu_SZ_Caudate_FU))
  d <- as.numeric(unlist(dfw$GABA_SZ_Caudate_FU))

  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]

  hc <- data.frame()

  dfw <- list(HC = data.frame(Glu_SZBL_Caudate = a1, GABA_SZBL_Caudate = b1),
              SZ = data.frame(Glu_SZ_Caudate = c1, GABA_SZ_Caudate = d1))

  ccCa <- cocor(~Glu_SZBL_Caudate + GABA_SZBL_Caudate | Glu_SZ_Caudate + GABA_SZ_Caudate, dfw, alternative = 'two.sided')


  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)

  a <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_BL))
  b <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_BL))
  c <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_FU))
  d <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_FU))

  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]

  hc <- data.frame()

  dfw <- list(HC = data.frame(Glu_SZBL_Thalamus = a1, GABA_SZBL_Thalamus = b1),
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))

  ccTh <- cocor(~Glu_SZBL_Thalamus + GABA_SZBL_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')


###################
### Right #########
###################
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_SZ_Caudate_BL_R))
  b <- as.numeric(unlist(dfw$GABA_SZ_Caudate_BL_R))
  c <- as.numeric(unlist(dfw$Glu_SZ_Caudate_FU_R))
  d <- as.numeric(unlist(dfw$GABA_SZ_Caudate_FU_R))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_SZBL_Caudate = a1, GABA_SZBL_Caudate = b1),
              SZ = data.frame(Glu_SZ_Caudate = c1, GABA_SZ_Caudate = d1))
  
  ccCa <- cocor(~Glu_SZBL_Caudate + GABA_SZBL_Caudate | Glu_SZ_Caudate + GABA_SZ_Caudate, dfw, alternative = 'two.sided')
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_BL_R))
  b <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_BL_R))
  c <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_FU_R))
  d <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_FU_R))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_SZBL_Thalamus = a1, GABA_SZBL_Thalamus = b1),
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))
  
  ccTh <- cocor(~Glu_SZBL_Thalamus + GABA_SZBL_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')
  
  ###################
  ### Left #########
  ###################
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_SZ_Caudate_BL_L))
  b <- as.numeric(unlist(dfw$GABA_SZ_Caudate_BL_L))
  c <- as.numeric(unlist(dfw$Glu_SZ_Caudate_FU_L))
  d <- as.numeric(unlist(dfw$GABA_SZ_Caudate_FU_L))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_SZBL_Caudate = a1, GABA_SZBL_Caudate = b1),
              SZ = data.frame(Glu_SZ_Caudate = c1, GABA_SZ_Caudate = d1))
  
  ccCa <- cocor(~Glu_SZBL_Caudate + GABA_SZBL_Caudate | Glu_SZ_Caudate + GABA_SZ_Caudate, dfw, alternative = 'two.sided')
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_BL_L))
  b <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_BL_L))
  c <- as.numeric(unlist(dfw$Glu_SZ_Thalamus_FU_L))
  d <- as.numeric(unlist(dfw$GABA_SZ_Thalamus_FU_L))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_SZBL_Thalamus = a1, GABA_SZBL_Thalamus = b1),
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))
  
  ccTh <- cocor(~Glu_SZBL_Thalamus + GABA_SZBL_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')
  
  
  ######################
  #### NAA GPC ratio ###
  ######################
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glu','GPC','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAA_HC_Caudate_BL))
  b <- as.numeric(unlist(dfw$GPC_HC_Caudate_BL))
  c <- as.numeric(unlist(dfw$NAA_SZ_Caudate_BL))
  d <- as.numeric(unlist(dfw$GPC_SZ_Caudate_BL))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(NAA_HC_Caudate = a1, GPC_HC_Caudate = b1), 
              SZ = data.frame(NAA_SZ_Caudate = c1, GPC_SZ_Caudate = d1))
  
  ccCa <- cocor(~NAA_HC_Caudate + GPC_HC_Caudate | NAA_SZ_Caudate + GPC_SZ_Caudate, dfw, alternative = 'two.sided')
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glu','GPC','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAA_HC_Thalamus_BL))
  b <- as.numeric(unlist(dfw$GPC_HC_Thalamus_BL))
  c <- as.numeric(unlist(dfw$NAA_SZ_Thalamus_BL))
  d <- as.numeric(unlist(dfw$GPC_SZ_Thalamus_BL))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Thalamus = a1, GABA_HC_Thalamus = b1), 
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))
  
  ccTh <- cocor(~Glu_HC_Thalamus + GABA_HC_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')
  
  ##################
  ### Right ########
  ##################
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAA_HC_Caudate_BL_R))
  b <- as.numeric(unlist(dfw$GPC_HC_Caudate_BL_R))
  c <- as.numeric(unlist(dfw$NAA_SZ_Caudate_BL_R))
  d <- as.numeric(unlist(dfw$GPC_SZ_Caudate_BL_R))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(NAA_HC_Caudate = a1, GPC_HC_Caudate = b1), 
              SZ = data.frame(NAA_SZ_Caudate = c1, GPC_SZ_Caudate = d1))
  
  ccCa <- cocor(~NAA_HC_Caudate + GPC_HC_Caudate | NAA_SZ_Caudate + GPC_SZ_Caudate, dfw, alternative = 'two.sided')
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAA_HC_Thalamus_BL_R))
  b <- as.numeric(unlist(dfw$GPC_HC_Thalamus_BL_R))
  c <- as.numeric(unlist(dfw$NAA_SZ_Thalamus_BL_R))
  d <- as.numeric(unlist(dfw$GPC_SZ_Thalamus_BL_R))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Thalamus = a1, GABA_HC_Thalamus = b1), 
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))
  
  ccTh <- cocor(~Glu_HC_Thalamus + GABA_HC_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')
  
  
  ##################
  ### Left ########
  ##################
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAA_HC_Caudate_BL_L))
  b <- as.numeric(unlist(dfw$GPC_HC_Caudate_BL_L))
  c <- as.numeric(unlist(dfw$NAA_SZ_Caudate_BL_L))
  d <- as.numeric(unlist(dfw$GPC_SZ_Caudate_BL_L))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(NAA_HC_Caudate = a1, GPC_HC_Caudate = b1), 
              SZ = data.frame(NAA_SZ_Caudate = c1, GPC_SZ_Caudate = d1))
  
  ccCa <- cocor(~NAA_HC_Caudate + GPC_HC_Caudate | NAA_SZ_Caudate + GPC_SZ_Caudate, dfw, alternative = 'two.sided')
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAA_HC_Thalamus_BL_L))
  b <- as.numeric(unlist(dfw$GPC_HC_Thalamus_BL_L))
  c <- as.numeric(unlist(dfw$NAA_SZ_Thalamus_BL_L))
  d <- as.numeric(unlist(dfw$GPC_SZ_Thalamus_BL_L))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Thalamus = a1, GABA_HC_Thalamus = b1), 
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))
  
  ccTh <- cocor(~Glu_HC_Thalamus + GABA_HC_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')
  
  #######################
  #### NAA NAAG ratio ###
  #######################
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAA_HC_Thalamus_BL))
  b <- as.numeric(unlist(dfw$NAAG_HC_Thalamus_BL))
  c <- as.numeric(unlist(dfw$NAA_SZ_Thalamus_BL))
  d <- as.numeric(unlist(dfw$NAAG_SZ_Thalamus_BL))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(Glu_HC_Thalamus = a1, GABA_HC_Thalamus = b1), 
              SZ = data.frame(Glu_SZ_Thalamus = c1, GABA_SZ_Thalamus = d1))
  
  ccTh <- cocor(~Glu_HC_Thalamus + GABA_HC_Thalamus | Glu_SZ_Thalamus + GABA_SZ_Thalamus, dfw, alternative = 'two.sided')
  
  #######################
  #### NAAG L vs R   ###
  #######################
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfw <- as.data.frame(dfw)
  
  a <- as.numeric(unlist(dfw$NAAG_HC_Thalamus_BL_R))
  b <- as.numeric(unlist(dfw$NAAG_HC_Thalamus_BL_L))
  c <- as.numeric(unlist(dfw$NAAG_SZ_Thalamus_BL_R))
  d <- as.numeric(unlist(dfw$NAAG_SZ_Thalamus_BL_L))
  
  a1 <- a[!is.na(a) & !is.na(b)]
  b1 <- b[!is.na(a) & !is.na(b)]
  c1 <- c[!is.na(c) & !is.na(d)]
  d1 <- d[!is.na(c) & !is.na(d)]
  
  hc <- data.frame()
  
  dfw <- list(HC = data.frame(NAAG_HC_Thalamus_BL_R = a1, NAAG_HC_Thalamus_BL_L = b1), 
              SZ = data.frame(NAAG_SZ_Thalamus_BL_R = c1, NAAG_SZ_Thalamus_BL_L = d1))
  
  ccTh <- cocor(~NAAG_HC_Thalamus_BL_R + NAAG_HC_Thalamus_BL_L | NAAG_SZ_Thalamus_BL_R + NAAG_SZ_Thalamus_BL_L, dfw, alternative = 'two.sided')
  
  
}

#############################
### Group Level R/NR.     ###
#############################

if (group_r_nr==T){
  
  df$group <- relevel(as.factor(df$group),'HC')
  df1 <- df %>% filter(!is.na(group))
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #ThGln <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  #ThGpc <- tidy(lm(data = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GPC ~ group*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')
  
  
  CaGlu <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #CaGln <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaNAAG <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC ~ group*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC.Cho ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAAG,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  
  Mall_R_Rem <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupR:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_Rem <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_R_NR <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupNR:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_NR <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupNR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  
  CaGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaGlu, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  Ca <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(Ca, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  
  fv <- fitted(Ca)
  rsq <- (residuals(Ca))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  # 
  M <- glm(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'),  mI ~ group*hemi + sex + scale(GMrat), family = inverse.gaussian(link = "1/mu^2"))
  emm <- emmeans(M, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  
  Th <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(Th, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  fv <- fitted(Th)
  rsq <- (residuals(Th))^2
  valid_indices <- which(fv > 0 & rsq > 0)
  ln_res_sq <- log(rsq[valid_indices])
  ln_fitted <- log(fv[valid_indices])
  park_model <- lm(ln_res_sq ~ ln_fitted)
  summary(park_model)
  
  
  CaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaGABA, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")

  CaGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaGluGln, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  ThGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(ThGABA, ~ group | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # CaGln <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))
  # emm <- emmeans(CaGln, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # # Convert to data frame
  # emm_df <- as.data.frame(emm)
  # pairs(emm,adjust = "fdr")
 
  CaGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaGlu, ~ group | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
   
  CaNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaNAA, ~ group | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # CaNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))
  # emm <- emmeans(CaNAA, ~ group | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # # Convert to data frame
  # emm_df <- as.data.frame(emm)
  # pairs(emm,adjust = "fdr")
  
  ThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(ThNAAG, ~ group | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
}

if (clinical==T){
  
  
  #add in sex from outside spreadsheet
  # 2026-06-02 AndyP changed to sheet 3 was reading incorrect sheet
  sds <- read_excel('7T.DARES.baselines.xlsx', sheet = 3) %>% group_by(RECID) %>% slice(1) %>% ungroup()
  sds <- sds %>% rename(id = RECID)
  sds$id <- as.character(sds$id)
  

  dfc <- left_join(df,sds,by='id')
  
  
  
  ############################
  #add in remitter status from outside spreadsheet
  sds2 <- read_excel('BPRS_items_MIKE.xlsx')
  sds2 <- sds2 %>% filter(Timepoint == 0 | Timepoint == 6) %>% mutate(timepoint = case_when(Timepoint == 0 ~ 'BL',Timepoint == 6 ~ 'FU'))
  sds2 <- sds2 %>% rename(id = RECID)
  sds2$id <- as.character(sds2$id)
  
  dfc <- left_join(dfc,sds2,by=c('id','timepoint'))
  
  sds3 <- read_excel('CGI values.xlsx') %>% mutate(timepoint = case_when(TP == 0 ~ 'BL',TP == 1 ~ 'FU'))
  sds3 <- sds3 %>% rename(id = RECID)
  sds3$id <- as.character(sds3$id)
  
  dfc <- left_join(dfc, sds3, by=c('id','timepoint'))

  dfc$timepoint <- relevel(factor(dfc$timepoint),ref = 'BL')
  dfc$group <- relevel(factor(dfc$Remitter_Status),ref = 'R')
  dfc$sex <- relevel(factor(dfc$sex.x),ref = "F")
  dfc$dup_mo <- dfc$`DUP (months)`
  
  dfc <- dfc %>% mutate(possx_sc = (POSSX.x), sev_sc = (Severity), imp_sc = (Improvement), 
                      wm_sc = scale(`Working Memory`), mccb_sc = scale(`MCCB Overall`), vl_sc = scale(`Verbal Learning`),
                      ps_sc = scale(`Processing Speed`),sct_sc = scale(`Strauss-Carpenter Total`),rpc_sc = scale(`Reasoning-Problem Solving`),
                      sc_sc = scale(`Social Cognition`))
  
  dfc <- dfc %>% dplyr::select(!`Working Memory` & !`Verbal Learning` & !`Strauss-Carpenter Total` & !`Processing Speed` & !`Reasoning-Problem Solving` & !`Social Cognition` & !`MCCB Overall`)
  dfrem <- dfc %>% dplyr::select(id,sex,timepoint,doi_m,wm_sc,possx_sc,sev_sc,imp_sc,mccb_sc,vl_sc,sc_sc,ps_sc,sct_sc,rpc_sc)
  dfrem <- dfrem %>% group_by(id,sex,timepoint) %>% slice(1) %>% ungroup()
  dfrem <- dfrem %>% pivot_wider(id_cols = c('id','sex','doi_m'), names_from = 'timepoint', values_from = c('possx_sc','sev_sc','imp_sc'))
  dfrem <- dfrem %>% mutate(dpossx = (possx_sc_FU - possx_sc_BL)/possx_sc_BL,
                            dsev = (sev_sc_FU - sev_sc_BL),
                            dimp = (imp_sc_FU - imp_sc_BL))
  ddfc <- dfc %>% pivot_wider(id_cols = c('id','roi','hemi','sex'), names_from = 'timepoint', values_from = c('GMrat','Glu','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'),values_fill = 0)
  ddfc <- ddfc %>% mutate(dGlu = Glu_FU - Glu_BL,
                      dGPC.Cho = GPC.Cho_FU - GPC.Cho_BL,
                      dGABA = GABA_FU - GABA_BL,
                      dNAA = NAA_FU - NAA_BL,
                      dmI = mI_FU - mI_BL,
                      dGln = Gln_FU - Gln_BL,
                      dNAAG = NAAG_FU - NAAG_BL,
                      dGlu.Gln = Glu.Gln_FU - Glu.Gln_BL)
  dfrem <- dfrem %>% dplyr::select(!sex)
  ddfc <- inner_join(ddfc,dfrem,by=c('id'))
  
  dfc$possx_sc <- scale(dfc$possx_sc)
  dfc$dup_mo_sc <- scale(Winsorize(dfc$dup_mo, val = quantile(dfc$dup_mo,probs = c(0.05,0.95),na.rm=TRUE)))
  
  #df1 <- df #%>% filter(doi_yr >= median(doi_yr,na.rm=TRUE))
  dfc <- dfc %>% ungroup() %>% mutate(doi_group = case_when(doi_m >= median(doi_m,na.rm=TRUE) ~ 'high',
                                              doi_m < median(doi_m,na.rm=TRUE) ~ 'low',
                                              group_level == 'HC' ~ 'HC'))
  dfc$doi_group <- relevel(as.factor(dfc$doi_group),ref='HC')
  
  dfc$curr_var <- dfc$imp_sc
  
  ThGlu <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Thalamus'), Glu ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Thalamus'), GABA ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Thalamus'), mI ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGluGln <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Thalamus'), Glu.Gln ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAA <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Thalamus'), NAA ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  ThNAAG <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Thalamus'), NAAG ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThGpc <- tidy(lm(data = dfc %>% filter(roi == 'Thalamus'), GPC ~ curr_var*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  Mth <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  #############################
  ####### Caudate #############
  #############################
  
  CaGlu <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), Glu ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), GABA ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), mI ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CaGluGln <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), Glu.Gln ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaNAAG <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), NAAG ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), NAA ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = dfc %>% filter(roi == 'Caudate'), GPC ~ curr_var*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), GPC.Cho ~ curr_var*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAAG,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')

  Mall <- rbind(Mth,Mca) %>% filter(term %in% c('curr_var')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_r <- rbind(Mth,Mca) %>% filter(term %in% c('curr_var:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)

  
  # examine emmeans
  M <- lmerTest::lmer(data = dfc %>% filter(roi == 'Caudate'), mI ~ curr_var*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ curr_var | hemi, data = dfc %>% filter(roi == 'Caudate'))
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  M <- lmerTest::lmer(data = dfc %>% filter(roi == 'Thalamus'), NAAG ~ curr_var*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ curr_var | hemi, data = dfc %>% filter(roi == 'Thalamus'))
  pairs(emm,adjust = "fdr")
  
  ddfc$id <- as.numeric(ddfc$id)
  ddfc <- ddfc %>% filter(id < 10000) # simple way to get rid of HC
  ddfc$id <- as.character(ddfc$id)
  
  ddfc$curr_var <- scale(ddfc$dsev)
  #ddfc <- ddfc %>% filter(dpossx < 0)
  
  ThGlu <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Thalamus'), dGlu ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Thalamus'), dGABA ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Thalamus'), dmI ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGluGln <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Thalamus'), dGlu.Gln ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAA <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Thalamus'), dNAA ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  ThNAAG <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Thalamus'), dNAAG ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  #ThGpc <- tidy(lm(data = ddfc %>% filter(roi == 'Thalamus'), dGPC ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU))) %>% mutate(metabolite = 'GPC')
  Mth <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')
  
  #############################
  ####### Caudate #############
  #############################
  
  CaGlu <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Caudate'), dGlu ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Caudate'), dGABA ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU)+ (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Caudate'), dmI ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CaGluGln <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Caudate'), dGlu.Gln ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaNAAG <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Caudate'), dNAAG ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Caudate'), dNAA ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  #CaGpc <- tidy(lm(data = ddfc %>% filter(roi == 'Caudate'), dGPC ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = ddfc %>% filter(roi == 'Caudate'), dGPC.Cho ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAAG,CaNAA,CaCho) %>% mutate(roi = 'Caudate')
  
  Mall <- rbind(Mth,Mca) %>% filter(term %in% c('curr_var')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  Mall_r <- rbind(Mth,Mca) %>% filter(term %in% c('curr_var:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% filter(pfdr < 0.05)
  
  
}

if (handedness_group == T){
  
  dfR <- df %>% filter(group_level == 'SZ' & hand %in% c('L','R'))
  dfR$hand <- relevel(as.factor(dfR$hand),'R')
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), Glu ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), GABA ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), mI ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGln <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), Gln ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), Glu.Gln ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), NAAG ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), NAA ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  ThGpc <- tidy(lm(data = dfR %>% filter(roi == 'Thalamus'), GPC ~ hemi*hand + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), GPC.Cho ~ group_level + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGln,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  
  CaGlu <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), Glu ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), GABA ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), mI ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CaGln <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), Gln ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), Glu.Gln ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaNAAG <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), NAAG ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), NAA ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = dfR %>% filter(roi == 'Caudate'), GPC ~ hemi*hand + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), GPC.Cho ~ hemi*hand + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGln,CaGluGln,CaNAAG,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  Mall_R_LH <- rbind(Mth1,Mca1) %>% filter(term %in% c('hemiR:handL')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_LH <- rbind(Mth1,Mca1) %>% filter(term %in% c('handL')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_R <- rbind(Mth1,Mca1) %>% filter(term %in% c('hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  #Mall_SZ <- rbind(Mth1,Mca1) %>% filter(term %in% c('group_levelSZ')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  
  
  library(emmeans)
  
  # examine emmeans
  MCaGlu <- lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), Glu ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGlu, ~ hemi, data = dfR %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGlu <- lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), Glu ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGlu, ~ hemi, data = dfR %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), Glu.Gln ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGluGln, ~ hemi, data = dfR %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), NAAG ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAAG, ~ hemi, data = dfR %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), NAAG ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAAG, ~ hemi, data = dfR %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAA <- lmerTest::lmer(data = dfR %>% filter(roi == 'Thalamus'), NAA ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAA, ~ hemi, data = dfR %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # # Create a customized plot
  # ggplot(emm_dfR, aes(x = group_level, y = emmean, color = hemi, group=hemi)) +
  #   geom_point(size = 3) +
  #   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  #   theme_minimal() +
  #   labs(title = "Model-Predicted Means by Group",
  #        y = "Estimated Marginal Mean",
  #        x = "Group")
  MCaGluGln <- lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), Glu.Gln ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGluGln, ~ hemi, data = dfR %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  MCaGABA <- lmerTest::lmer(data = dfR %>% filter(roi == 'Caudate'), GABA ~ hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGABA, ~ hemi, data = dfR %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_dfR <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
}

if (Figure_2){
  
  secondary_mets <- c('GPC','NAAG','GPC.Cho','mI','NAA')
  primary_mets <- c('Glu','GABA','Glu.Gln')
  
  
  dodge_width = 0.8
  df <- df %>% select(id,group_level,GMrat,GPC,Glu,GPC.Cho,GABA,NAA,mI,Gln,NAAG,Glu.Gln,roi,hemi,timepoint)
  dfL <- df %>% filter(timepoint == 'BL') %>% pivot_longer(cols = c('GMrat','GPC','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfL <- dfL %>% filter(name != 'Gln' & name != 'GMrat')
  dfL$Group <- dfL$group_level
  
  ####### Left Caudate Primary Metabolites ########
  pdf('Figure_2A_L_Caudate_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% secondary_mets) & hemi == 'L'), aes(x = name, y = value)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "***",              # Custom text or star
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  ####### Right Caudate Primary Metabolites ########
  pdf('Figure_2B_R_Caudate_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% secondary_mets) & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "**",              # Custom text or star
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "*",
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "***",
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  ####### Left Thalamus Primary Metabolites ########
  pdf('Figure_2C_L_Thalamus_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% secondary_mets) & hemi == 'L'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "*",              # Custom text or star
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  ####### Right Thalamus Primary Metabolites ########
  pdf('Figure_2D_R_Thalamus_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% secondary_mets) & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",              # Custom text or star
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "**",
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  
  ####### Left Caudate Secondary Metabolites ########
  pdf('Figure_2E_L_Caudate_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% primary_mets) & !name == 'NAAG' & hemi == 'L'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",              # Custom text or star
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "**",
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "p=0.066",
      vjust = -0.75,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 4 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 4 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "*",
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  ####### Right Caudate Secondary Metabolites ########
  pdf('Figure_2F_R_Caudate_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% primary_mets) & !name == 'NAAG' & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "*",              # Custom text or star
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 4 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 4 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  ####### Left Thalamus Secondary Metabolites ########
  pdf('Figure_2G_L_Thalamus_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% primary_mets) & name != 'GPC' & hemi == 'L'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "*",              # Custom text or star
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "***",
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 6,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 4 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 4 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "***",
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  ####### Right Thalamus Secondary Metabolites ########
  pdf('Figure_2H_R_Thalamus_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% primary_mets) & name != 'GPC' & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 1 - (dodge_width / 4),  # Centers on the left bar (approx 0.8)
      xmax = 1 + (dodge_width / 4),  # Centers on the right bar (approx 1.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",              # Custom text or star
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 2 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 2 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 3 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 3 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "N.S.",
      vjust = -0.5,
      tip_length = 0.03
    ) +
    
    # Bracket 2: Compare Control vs Treat INSIDE "Site B" (X = 2)
    geom_signif(
      color = "black",
      textsize = 4,          # Optional: Adjusts text size of the labels
      linewidth = 1,       # Optional: Adjusts the thickness of the bracket lines
      xmin = 4 - (dodge_width / 4),  # Centers on Site B's left bar (approx 1.8)
      xmax = 4 + (dodge_width / 4),  # Centers on Site B's right bar (approx 2.2)
      y_position = 2.5,               # Height of the bracket
      annotation = "*",
      tip_length = 0.03
    ) + theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
}

if (Creatine_Check){
  
  
  dg <- read_csv('13MP20200207_LCMv2fixidx_Raw.csv')
  dg <- dg %>% separate_wider_delim(cols = ld8,delim="_",names=c("id","dateNumeric"),cols_remove=FALSE)
  dg$dateNumeric <- as.numeric(dg$dateNumeric)
  dg <- dg %>% group_by(id,visitnum,label) %>% slice(1) %>% ungroup()
  dg <- dg %>% mutate(Cr = case_when(Cr.SD > 20 | Cr < 0.01 | Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ Cr))
  dg <- dg %>% filter(age >= 18) %>%
    filter(label %in% c('L Thalamus','R Thalamus','L Caudate','R Caudate')) %>% 
                        group_by(id,label) %>% mutate(nsess = 1:n()) %>%
             filter(nsess == 1) %>% mutate(timepoint = 'BL')
  
  # Match ROI names, create date value
  dg <- dg %>% 
    separate(label, c('hemi','roi'), ' ', convert=TRUE) %>%
    mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
    mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
    mutate(hemi = ifelse(hemi=='L', 'L', ifelse(hemi=='R', 'R', NA))) %>% # Fixed here
    mutate(region = paste0(hemi,' ', roi))
  dg <- dg %>% select(id,roi,hemi,timepoint,Cr, Cr.SD)
                      
  # hc_mike <- read_excel('13MP20200207_LCMv2fixidx_Mike.xlsx') %>% separate_wider_delim(cols = RECID, delim = "_", names = c('id','date'))
  # #rm(met_out1)
  # Load SZ MRSI data
  szmet_orig <- readxl::read_xlsx('sarpal_mrsi_original_07062026.xlsx', sheet = 1) %>%
    rename(scan_date = Scan_date) %>% mutate(source = 'orig')
  #szmet_orig <- szmet_orig %>% group_by(RECID) %>% separate_wider_delim(cols = region, delim = " ",names = c('hemisphere','roi'), cols_remove = TRUE,too_few = "align_start") %>% ungroup()
  szmet_orig <- szmet_orig %>% rename(roi = 'region')
  # Load SZ MRSI data
  szmet_orig <- readxl::read_xlsx('sarpal_mrsi_original_07062026.xlsx', sheet = 1) %>%
    rename(scan_date = Scan_date) %>% mutate(source = 'orig')
  #szmet_orig <- szmet_orig %>% group_by(RECID) %>% separate_wider_delim(cols = region, delim = " ",names = c('hemisphere','roi'), cols_remove = TRUE,too_few = "align_start") %>% ungroup()
  szmet_orig <- szmet_orig %>% rename(roi = 'region') %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                                                                     timepoint == 2 ~ 'FU'))


  # Match ROI names, create date value
  szmet_new <- szmet_orig %>% 
    separate(roi, c('hemi','roi'), ' ', convert=TRUE) %>%
    mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
    mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
    mutate(dateNumeric = as.numeric(as.POSIXct(scan_date, format="%Y-%m-%d")),
         hemi = ifelse(hemi=='left', 'L', ifelse(hemi=='right', 'R', NA))) %>%
         mutate(region = paste0(hemi,' ', roi))
  
  # Match ROI names, create date value
  szmet_new <- szmet_orig %>% 
    separate(roi, c('hemi','roi'), ' ', convert=TRUE) %>%
    mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
    mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
    mutate(dateNumeric = as.numeric(as.POSIXct(scan_date, format="%Y-%m-%d")),
           hemi = ifelse(hemi=='left', 'L', ifelse(hemi=='right', 'R', NA))) %>%
    mutate(region = paste0(hemi,' ', roi))
  
  szmet_new <- szmet_new %>% select(RECID,timepoint,roi,hemi,Cre,`Cre %SD`)
  szmet_new <- szmet_new %>% rename(id = RECID, Cr = Cre, Cr.SD = `Cre %SD`)
  szmet_new$id <- as.character(szmet_new$id)
  szmet_new <- szmet_new %>% mutate(Cr = case_when(Cr.SD > 20 | Cr < 0.01 | Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ Cr))
  
  df <- df %>% filter((group_level == 'HC' & timepoint == 'BL') | group_level == 'SZ')
  
  
  dg <- rbind(dg,szmet_new)
  df <- left_join(df,dg,by=c('id','timepoint','roi','hemi'))
  
  df <- df %>% mutate(Cr_gamadj1 = case_when(Cr.SD > 20 | Cr < 0.01 | Cr > 5*sd(Cr,na.rm=T) ~ NA_real_, TRUE ~ Cr))
  df$Cr_gamadj1 <- scale(Winsorize(df$Cr_gamadj, val = quantile(df$Cr_gamadj,probs = c(0.05,0.95),na.rm=TRUE)))
  
  df <- df %>% group_by(id,roi,hemi,timepoint) %>% slice(1) %>% ungroup()
  
  MCrCa <- lmerTest::lmer(data= df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Cr_gamadj1 ~ group_level*hemi + scale(GMrat) + sex + (1|id))
  MCrTh <- lmerTest::lmer(data= df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Cr_gamadj1 ~ group_level*hemi + scale(GMrat) + sex + (1|id))
  
  MCrCaG <- lmerTest::lmer(data= df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Cr_gamadj1 ~ group*hemi + scale(GMrat) + sex + (1|id))
  MCrThG <- lmerTest::lmer(data= df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Cr_gamadj1 ~ group*hemi + scale(GMrat) + sex + (1|id))
  
  
  MCrCaL <- lmerTest::lmer(data= df %>% filter(roi == 'Caudate'), Cr_gamadj1 ~ condition*hemi + scale(GMrat) + sex + (1|id))
  MCrThL <- lmerTest::lmer(data= df %>% filter(roi == 'Thalamus'), Cr_gamadj1 ~ condition*hemi + scale(GMrat) + sex + (1|id))
  
  emm <- emmeans(MCrCa, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  emm <- emmeans(MCrTh, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  emm <- emmeans(MCrCaL, ~ condition | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  emm <- emmeans(MCrThL, ~ condition | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  emm <- emmeans(MCrCaG, ~ group | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  emm <- emmeans(MCrThG, ~ group | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ group_level*hemi + scale(GMrat) + Cr_gamadj*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ group_level*hemi + scale(GMrat) + Cr_gamadj*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ group_level*hemi + scale(GMrat) + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Gln ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + scale(GMrat) + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group_level*hemi + scale(GMrat) + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAA ~ group_level*hemi + scale(GMrat) + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  #ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GPC ~ group_level*hemi + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group_level + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')
  
  
  # CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  # CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  # CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  # #CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  # CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  # #CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  # CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # # convergence issues in this lmer model, low variance at the subject level so just use lm
  # CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC ~ group_level*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  # CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC.Cho ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  # Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  # 
  Mall_R_SZ <- rbind(Mth1) %>% filter(term %in% c('group_levelSZ:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr) %>% filter(pfdr < 0.05)
  Mall_SZ <- rbind(Mth1) %>% filter(term %in% c('group_levelSZ')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr) %>% filter(pfdr < 0.05)
  
  
  # examine emmeans
  MCaGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'),  Glu ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCaGlu, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCaGPC.Cho <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'),  GPC.Cho ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCaGPC.Cho, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCaGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThGlu, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThGluGln, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThGluGln, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThGluGln, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThmI <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThmI, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCamI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCamI, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans, yes there is more in L in HC and more in R in SZ
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThNAAG, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThNAA, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCaNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCaNAA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCaNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCaNAA, ~ hemi | group_level, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCamI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCamI, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  # # Create a customized plot
  # ggplot(emm_df, aes(x = group_level, y = emmean, color = hemi, group=hemi)) +
  #   geom_point(size = 3) +
  #   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  #   theme_minimal() +
  #   labs(title = "Model-Predicted Means by Group",
  #        y = "Estimated Marginal Mean",
  #        x = "Group")
  MCaGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCaGluGln, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  MCaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MCaGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  df$Group <- df$group_level
  pdf('Cr_baseline_Comparison.pdf',height=5,width=6)
  gg1 <- ggplot(data = df, aes(x=Group,y = Cr_gamadj,color = Group)) + geom_jitter() + 
    geom_boxplot(color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    theme_minimal() + facet_grid(roi~hemi) + ylim(c(0,1000)) +    
    theme(axis.text.y = element_text(size = 14),
                                axis.text.x = element_text(size = 14),
                                legend.text = element_text(size = 14),
                                legend.title = element_text(size = 16),
                                axis.title.x = element_text(size = 16),
                                axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  df$Group <- df$timepoint
  pdf('Cr_baseline_Comparison_longitudinal.pdf',height=5,width=6)
  gg1 <- ggplot(data = df %>% filter(group_level == 'SZ'), aes(x=Group,y = Cr_gamadj)) + geom_jitter() + 
    geom_boxplot(color = "black",outlier.shape = NA,notch = T,linewidth = 0.75,alpha = 0.5,fatten = 1) +
    theme_minimal() + facet_grid(roi~hemi) + ylim(c(0,1000)) +   
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + ylab('Cr') + xlab('Timepoint')
  print(gg1)
  dev.off()
  
  df$Group <- df$group
  pdf('Cr_baseline_Comparison_R_NR.pdf',height=5,width=6)
  gg1 <- ggplot(data = df %>% filter(group_level == 'SZ' & !is.na(group) & timepoint == 'BL'), aes(x=Group,y = Cr_gamadj)) + geom_jitter() + 
    geom_boxplot(color = "black",outlier.shape = NA,notch = F,linewidth = 0.75,alpha = 0.5,fatten = 1) +
    theme_minimal() + facet_grid(roi~hemi) + ylim(c(0,1000)) +  
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + ylab('Cr') + xlab('Timepoint')
  print(gg1)
  dev.off()
  
  
}

if (loglink_GLM==T){
  
  # df1 <- df %>% mutate(Glu = Glu*Cr_gamadj,
  #                     GABA = GABA*Cr_gamadj,
  #                     mI = mI*Cr_gamadj,
  #                     Glu.Gln = Glu.Gln*Cr_gamadj,
  #                     NAAG = NAAG*Cr_gamadj,
  #                     NAA = NAA*Cr_gamadj)
  # df1 <- df1 %>% ungroup()
  # df1 = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL')
  # # 2026-06-04 Will need to add hemi eventually
  # ThGlu <- tidy(glm(data = df1, Glu ~ group_level*hemi +  offset(log(Cr_gamadj)), family = gaussian(link = "log"))) %>% mutate(metabolite = 'Glu')
  # ThGABA <- tidy(glm(data = df1, GABA ~ group_level*hemi+  offset(log(Cr_gamadj)), family = gaussian(link = "log"))) %>% mutate(metabolite = 'GABA')
  # ThGluGln <- tidy(glm(data = df1, Glu.Gln ~ group_level*hemi+  offset(log(Cr_gamadj)), family = gaussian(link = "log"))) %>% mutate(metabolite = 'GluGln')
  # ThNAAG <- tidy(glm(data = df1, NAAG ~ group_level*hemi +  offset(log(Cr_gamadj)), family = gaussian(link = "log"))) %>% mutate(metabolite = 'NAAG')
  # ThNAA <- tidy(glm(data = df1, NAA ~ group_level*hemi + offset(log(Cr_gamadj)), family = gaussian(link = "log"))) %>% mutate(metabolite = 'NAA')
  # Mth1 <- rbind(ThGlu,ThGABA,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')
  # 
  # Mall_R_SZ <- rbind(Mth1) %>% filter(term %in% c('group_levelSZ:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr)
  # Mall_SZ <- rbind(Mth1) %>% filter(term %in% c('group_levelSZ')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr)
  # 
  # 
  # # examine emmeans
  # MThNAAG <- glm(data = df1, NAAG ~ group_level*hemi +  offset(log(Cr_gamadj)), family = gaussian(link = "log"))
  # emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df1)
  # # Convert to data frame
  # emm_df <- as.data.frame(emm)
  # pairs(emm,adjust = "fdr")
  
  df1 <- df %>% mutate(Glu = Glu/GPC.Cho,
                       GABA = GABA/GPC.Cho,
                       mI = mI/GPC.Cho,
                       Glu.Gln = Glu.Gln/GPC.Cho,
                       NAAG = NAAG/GPC.Cho,
                       NAA = NAA/GPC.Cho)
  df1 <- df1 %>% ungroup()
  df1 = df1 %>% filter(roi == 'Thalamus' & timepoint == 'BL')
  # 2026-06-04 Will need to add hemi eventually
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df1, Glu ~ group_level*hemi + scale(GMrat) + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df1, GABA ~ group_level*hemi + scale(GMrat) + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df1, mI ~ group_level*hemi + scale(GMrat) + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGluGln <- tidy(lmerTest::lmer(data = df1, Glu.Gln ~ group_level*hemi + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df1, NAAG ~ group_level*hemi+ Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df1, NAA ~ group_level*hemi + Cr_gamadj1*group_level + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')

  Mall_R_SZ <- rbind(Mth1) %>% filter(term %in% c('group_levelSZ:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr)
  Mall_SZ <- rbind(Mth1) %>% filter(term %in% c('group_levelSZ')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr)
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = df1, NAAG ~ group_level*hemi+ Cr_gamadj1*group_level + (1|id))
  emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df1)
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
}

if (Figure_3){
  
  df1 <- df %>% filter(group_level == 'SZ' | (group_level == 'HC' & timepoint == 'BL'))
  
  secondary_mets <- c('GPC','NAAG','GPC.Cho','mI','NAA')
  primary_mets <- c('Glu','GABA','Glu.Gln')
  
  
  dodge_width = 0.8
  df1 <- df1 %>% select(id,group,GMrat,GPC,Glu,GPC.Cho,GABA,NAA,mI,Gln,NAAG,Glu.Gln,roi,hemi,timepoint) %>% filter(timepoint == 'BL')
  dfL <- df1 %>% pivot_longer(cols = c('GMrat','GPC','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfL <- dfL %>% filter(name != 'Gln' & name != 'GMrat' & !is.na(group))
  dfL$Group <- dfL$group
  
  pdf('Figure_3A_L_Caudate_BL.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% secondary_mets) & hemi == 'L'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Figure_3A_R_Caudate_BL.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% secondary_mets) & hemi == 'R'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Figure_3A_L_Thalamus_BL.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% secondary_mets) & hemi == 'L'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Figure_3A_R_Thalamus_BL.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% secondary_mets) & hemi == 'R'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Figure_3A_L_Caudate_BL_Secondary.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% primary_mets) & !name %in% c('NAAG','mI') & hemi == 'L'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Figure_3A_R_Caudate_BL_Secondary.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% primary_mets) & !name %in% c('NAAG','mI','GPC') & hemi == 'R'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Figure_3A_L_Thalamus_BL_Secondary.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% primary_mets) & !name %in% c('GPC','mI') & hemi == 'L'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Figure_3A_R_Thalamus_BL_Secondary.pdf',height = 4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% primary_mets) & !name %in% c('GPC','NAA','GPC.Cho','mI') & hemi == 'R'), aes(x = name, y = value,color=Group)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
    geom_boxplot(aes(Group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + theme_minimal() + xlab('') + ylab('') +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  
  df <- df %>% mutate(xvar = case_when(group_level == 'HC' ~ 'HC',
                                       group_level == 'SZ' & timepoint == 'BL' ~ 'Baseline',
                                       group_level == 'SZ' & timepoint == 'FU' ~ 'Follow-Up'))
  
  df$xvar <- relevel(as.factor(df$xvar),'HC')
  
  df <- df[order(df$id, df$xvar), ]
  # position = position_jitter(width = 0.3, height = 0, seed=45)
  # ggplot(df %>% filter(!is.na(group) & hemi == 'R' & roi == 'Caudate'), aes(x=xvar,y=Glu,color=group)) + 
  #   #geom_line(data = df %>% filter(group_level == 'SZ'),aes(group = id),position = position, alpha = 0.5) +
  #   geom_point(position = position, size =0.8) +
  #   geom_line(data = df %>% filter(xvar %in% c('Baseline','Follow-Up') & !is.na(group) & hemi == 'R' & roi == 'Caudate'),position = position,aes(group = id),alpha = 0.5)
  # 
  pdf('Figure_3B_R_Thalamus_Glutamate_Longitudinal.pdf',height = 4, width = 5)
  # 1. Create a consistent random shift for each ID
  set.seed(42)
  id_shifts <- df %>% filter(roi == 'Thalamus' & hemi == 'R' & !is.na(group))  %>% 
    distinct(id) %>%
    mutate(x_shift = runif(n(), min = -0.07, max = 0.07)) # Sightly narrower jitter to clear room for boxes
  
  # 2. Build a custom x-axis layout
  df_plot <- df %>% filter(roi == 'Thalamus' & hemi == 'R') %>%
    filter(!is.na(group)) %>%
    left_join(id_shifts, by = "id") %>%
    mutate(
      # Map the x-axis positions explicitly
      x_num = case_when(
        group == 'HC' ~ 1,                                # HC goes on the far left
        group != 'HC' & xvar == 'Baseline' ~ 2,           # SZ/Clinical Baseline in the middle
        group != 'HC' & xvar == 'Follow-Up' ~ 3,          # SZ/Clinical Follow-Up on the right
        TRUE ~ NA_real_
      ),
      group_dodge = ifelse(group == 'HC', -0.18, 0.18),
      x_center = x_num + 0.12,
      x_jittered = x_num -0.12 + x_shift
    ) %>%
    # Filter out any rows that didn't get mapped (just to be safe)
    filter(!is.na(x_num))
  
  # 3. Plot with the boxplot overlay
  gg1 <- ggplot(df_plot, aes(y = Glu, color = group, fill = group)) + 
    # Layer 1: Spaghetti lines (drawn in the background via alpha)
    geom_line(
      data = df_plot %>% filter(xvar %in% c('Baseline', 'Follow-Up')),
      aes(x = x_jittered, group = id),
      alpha = 0.35,
      show.legend = FALSE
    ) +
    # Layer 2: Raw data points
    geom_point(aes(x = x_jittered), size = 0.8, alpha = 0.7) +
    # Layer 3: Summary Boxplots positioned exactly over the cluster centers
    geom_boxplot(
      aes(x = x_center, group = interaction(xvar, group)),
      width = 0.15,          # Controls how wide the boxes are
      alpha = 0.4,           # Makes the box fill semi-transparent so lines show through
      outlier.shape = NA,     # CRUCIAL: Hides duplicate outlier points since geom_point already shows them
      notch = F
    ) +
    scale_x_continuous(
      breaks = c(1, 2), 
      labels = c('Baseline', 'Follow-Up')
    ) +
    theme_minimal() + scale_color_manual(values = c("HC" = "#E41A1C", "NR" = "#4DAF4A", "R" = "#377EB8")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "NR" = "#4DAF4A", "R" = "#377EB8")) +
    theme_minimal()
  print(gg1)
  dev.off()
  
  
  pdf('Figure_3B_R_Caudate_NAA_Longitudinal.pdf',height = 4, width = 5)
  # 1. Create a consistent random shift for each ID
  set.seed(42)
  id_shifts <- df %>% filter(roi == 'Caudate' & hemi == 'R' & !is.na(group))  %>% 
    distinct(id) %>%
    mutate(x_shift = runif(n(), min = -0.07, max = 0.07)) # Sightly narrower jitter to clear room for boxes
  
  # 2. Build a custom x-axis layout
  df_plot <- df %>% filter(roi == 'Caudate' & hemi == 'R') %>%
    filter(!is.na(group)) %>%
    left_join(id_shifts, by = "id") %>%
    mutate(
      # Map the x-axis positions explicitly
      x_num = case_when(
        group == 'HC' ~ 1,                                # HC goes on the far left
        group != 'HC' & xvar == 'Baseline' ~ 2,           # SZ/Clinical Baseline in the middle
        group != 'HC' & xvar == 'Follow-Up' ~ 3,          # SZ/Clinical Follow-Up on the right
        TRUE ~ NA_real_
      ),
      group_dodge = ifelse(group == 'HC', -0.18, 0.18),
      x_center = x_num + 0.12,
      x_jittered = x_num -0.12 + x_shift
    ) %>%
    # Filter out any rows that didn't get mapped (just to be safe)
    filter(!is.na(x_num))
  
  # 3. Plot with the boxplot overlay
  gg1 <- ggplot(df_plot, aes(y = NAA, color = group, fill = group)) + 
    # Layer 1: Spaghetti lines (drawn in the background via alpha)
    geom_line(
      data = df_plot %>% filter(xvar %in% c('Baseline', 'Follow-Up')),
      aes(x = x_jittered, group = id),
      alpha = 0.35,
      show.legend = FALSE
    ) +
    # Layer 2: Raw data points
    geom_point(aes(x = x_jittered), size = 0.8, alpha = 0.7) +
    # Layer 3: Summary Boxplots positioned exactly over the cluster centers
    geom_boxplot(
      aes(x = x_center, group = interaction(xvar, group)),
      width = 0.15,          # Controls how wide the boxes are
      alpha = 0.4,           # Makes the box fill semi-transparent so lines show through
      outlier.shape = NA,     # CRUCIAL: Hides duplicate outlier points since geom_point already shows them
      notch = F
    ) +
    scale_x_continuous(
      breaks = c(1, 2), 
      labels = c('Baseline', 'Follow-Up')
    ) +
    theme_minimal() + scale_color_manual(values = c("HC" = "#E41A1C", "NR" = "#4DAF4A", "R" = "#377EB8")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "NR" = "#4DAF4A", "R" = "#377EB8")) +
    theme_minimal()
  print(gg1)
  dev.off()
  
}

if (Figure_4){
  secondary_mets <- c('GPC','NAAG','GPC.Cho','mI','NAA')
  primary_mets <- c('Glu','GABA','Glu.Gln')

  
  dodge_width = 0.8
  df <- df %>% select(id,group_level,GMrat,GPC,Glu,GPC.Cho,GABA,NAA,mI,Gln,NAAG,Glu.Gln,roi,hemi,timepoint,doigr)
  dfL <- df %>% filter(timepoint == 'BL') %>% pivot_longer(cols = c('GMrat','GPC','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
  dfL <- dfL %>% filter(name != 'Gln' & name != 'GMrat')
  dfL$Group <- dfL$doigr
  
  dfL$Group <- factor(dfL$Group, levels = c('HC','low','high'))
  
  ####### Left Caudate Primary Metabolites ########
  pdf('Figure_4A_L_Caudate_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% secondary_mets) & hemi == 'L'), aes(x = name, y = value)) + 
    geom_jitter(aes(color = Group),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
  
  ####### Right Caudate Primary Metabolites ########
  pdf('Figure_4A_R_Caudate_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% secondary_mets) & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
  
  ####### Left Thalamus Primary Metabolites ########
  pdf('Figure_4A_L_Thalamus_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% secondary_mets) & hemi == 'L'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
  
  ####### Right Thalamus Primary Metabolites ########
  pdf('Figure_4A_R_Thalamus_Primary.pdf',height=4, width = 5)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% secondary_mets) & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
  
  
  ####### Left Caudate Secondary Metabolites ########
  pdf('Figure_4A_L_Caudate_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% primary_mets) & !name == 'NAAG' & hemi == 'L'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
  
  ####### Right Caudate Secondary Metabolites ########
  pdf('Figure_4A_R_Caudate_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Caudate' & !(name %in% primary_mets) & !name == 'NAAG' & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
  
  ####### Left Thalamus Secondary Metabolites ########
  pdf('Figure_4A_L_Thalamus_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% primary_mets) & name != 'GPC' & hemi == 'L'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
  
  ####### Right Thalamus Secondary Metabolites ########
  pdf('Figure_4A_R_Thalamus_Secondary.pdf',height=4, width = 7)
  gg1 <- ggplot(dfL %>% filter(roi == 'Thalamus' & !(name %in% primary_mets) & name != 'GPC' & hemi == 'R'), aes(x = name, y = value, color = Group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
                size = 0.8, alpha = 0.8) +
    geom_boxplot(aes(group = interaction(name,Group)),color = "black",outlier.shape = NA,notch = T,linewidth = 0.75, fill = NA,fatten = 1) +
    ylim(c(0,2.75)) + xlab('metabolite') + ylab('concentration (A.U.)') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) + scale_color_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A")) +
    scale_fill_manual(values = c("HC" = "#E41A1C", "low" = "#003366", "high" = "#EBCC2A"))
  print(gg1)
  dev.off()
}

if (corr_heatmap ==T){
  
  library(ggcorrplot)
  
  dfc <- df %>% filter(timepoint == 'BL') %>% group_by(id,group_level,roi) %>%
    summarize(across(where(is.numeric), \(x) mean(x, na.rm=TRUE))) %>%
    ungroup()
  
  dfc_hc_ca <- dfc %>% filter(group_level == 'HC' & roi == 'Caudate')
  dfc_dz_ca <- dfc %>% filter(group_level == 'SZ' & roi == 'Caudate')
  dfc_hc_th <- dfc %>% filter(group_level == 'HC' & roi == 'Thalamus')
  dfc_dz_th <- dfc %>% filter(group_level == 'SZ' & roi == 'Thalamus')
  
  col2rem <- c('id','group_level','roi','hemi','doi_m','Cr_gamadj','age','POSSX','GMrat','nsess','age_sc','Gln')
  
  remNAAGca <- c('NAAG')
  remGPCth <- c('GPC')
  dfc_hc_ca <- dfc_hc_ca %>% select(-any_of(col2rem)) %>% select(-any_of(remNAAGca))
  dfc_dz_ca <- dfc_dz_ca %>% select(-any_of(col2rem)) %>% select(-any_of(remNAAGca))
  dfc_hc_th <- dfc_hc_th %>% select(-any_of(col2rem)) %>% select(-any_of(remGPCth))
  dfc_dz_th <- dfc_dz_th %>% select(-any_of(col2rem)) %>% select(-any_of(remGPCth))
  
  cm_hc_ca <- cor(dfc_hc_ca, use = "complete.obs")
  cm_sz_ca <- cor(dfc_dz_ca, use = "complete.obs")
  cm_hc_th <- cor(dfc_hc_th, use = "complete.obs")
  cm_sz_th <- cor(dfc_dz_th, use = "complete.obs")
  
  p_matrix_cm_hc_ca <- cor_pmat(cm_hc_ca, use = "pairwise.complete.obs")
  p_matrix_cm_sz_ca <- cor_pmat(cm_sz_ca, use = "pairwise.complete.obs")
  p_matrix_cm_hc_th <- cor_pmat(cm_hc_th, use = "pairwise.complete.obs")
  p_matrix_cm_sz_th <- cor_pmat(cm_sz_th, use = "pairwise.complete.obs")
  
  pdf('Figure_SX_Correlation_Matrix_HC_Caudate.pdf',height=5,width=5)
  gg1 <- ggcorrplot(cm_hc_ca, lab = TRUE, hc.order = TRUE, 
             hc.method = "complete",type = "lower",
             p.mat = p_matrix_cm_hc_ca,
             sig.level = 0.05,
             insig = "blank",
             lab_size = 4, title = 'Bilateral Caudate HC') + theme_minimal()
  gg1 <- gg1 + theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
  print(gg1)
  dev.off()
  
  pdf('Figure_SX_Correlation_Matrix_SZ_Caudate.pdf',height=5,width=5)  
  gg1 <- ggcorrplot(cm_sz_ca, lab = TRUE, hc.order = TRUE, 
             hc.method = "complete",type = "lower",
             p.mat = p_matrix_cm_sz_ca,
             sig.level = 0.05,
             insig = "blank",
             lab_size = 4, title = 'Bilateral Caudate SSD') + theme_minimal()
  gg1 <- gg1 + theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
  print(gg1)
  dev.off()
  pdf('Figure_SX_Correlation_Matrix_HC_Thalamus.pdf',height=5,width=5) 
  gg1 <- ggcorrplot(cm_hc_th, lab = TRUE, hc.order = TRUE, 
             hc.method = "complete",type = "lower",
             p.mat = p_matrix_cm_hc_th,
             sig.level = 0.05,
             insig = "blank",
             lab_size = 4, title = 'Bilateral Thalamus HC') + theme_minimal()
  gg1 <- gg1 + theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
  print(gg1)
  dev.off()
  pdf('Figure_SX_Correlation_Matrix_SZ_Thalamus.pdf',height=5,width=5) 
  gg1 <- ggcorrplot(cm_sz_th, lab = TRUE, hc.order = TRUE, 
             hc.method = "complete",type = "lower",
             p.mat = p_matrix_cm_sz_th,
             sig.level = 0.05,
             insig = "blank",
             lab_size = 4, title = 'Bilateral Thalamus SSD') + theme_minimal()
  gg1 <- gg1 + theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
  print(gg1)
  dev.off()
 
  df$Group <- df$group_level
  pdf('Supplemental_XA_Glu_vs_GABA_Thalamus.pdf',height = 4,width=6)
  gg1 <- ggplot(data = df %>% filter(timepoint == 'BL' & roi == 'Thalamus'), aes(x=GABA, y=Glu,color= Group,group=Group)) + 
    geom_point(size=1) + geom_smooth(method = 'lm') + theme_minimal() +     
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
  pdf('Supplemental_XA_GPC_vs_NAA_Thalamus.pdf',height = 4,width=6)
  gg1 <- ggplot(data = df %>% filter(timepoint == 'BL' & roi == 'Caudate'), aes(x=GPC, y=NAA,color= Group,group=Group)) + 
    geom_point(size=1) + geom_smooth(method = 'lm') + theme_minimal() +     
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  print(gg1)
  dev.off()
  
}

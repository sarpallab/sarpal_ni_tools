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

reload_new = F # 2026-07-02 post redoing GM
# same as reload, but loads in non gam-adjusted metabolites, can merge with reload df manually
# also does group analyses with non-gam-adjusted metabolites
# old, part of redoing GM for revision using uniform method
reload = T
reload_check_GM = F
longitudinal = F
group = T
hilowdoi = F
eibalance = F
group_r_nr = F
# will reload just Sarpal / CZ data
clinical = F
handedness_group = F

#The manuscript states that FDR correction was performed by accounting for metabolites within each ROI. 
# However, the statistical inference and biological interpretation are made across three ROIs. 
# Please state exactly which p-values were included in each FDR correction family. 
# In addition, please provide a more conservative sensitivity analysis in which FDR correction 
# includes all ROI-by-metabolite tests within each type of analysis, and report whether the main findings remain significant.

# macbook
basedir <- ('/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper')
# macmini
#basedir <- setwd('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper')
setwd(basedir)

if (reload_new == T){
  hc <- read_csv(paste0(basedir,'/HC/','metabolites_gamadj_long_HC_07012026.csv')) %>% mutate(group = 'HC', timepoint = 'BL',`DUP (months)` = NA)
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
  # df <- df %>% filter(region != 'right caudate' & group == 'HC') #%>% mutate(RECID = gsub("^(.*?)_.*", "\\1", RECID))
  # 
  # testRC <- read_xlsx('Copy of R_Caudate_Glu.xlsx') %>% rename(RECID = id)#%>% 
  #   #mutate(RECID = gsub("^(.*?)_.*", "\\1", RECID))
  
  #df <- inner_join(testRC,df, by=c('RECID'))
  
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
  df <- df %>% dplyr::select(RECID,doi_m,region,timepoint,age,sex,group,Scan_date,POSSX,Remitter_Status,GMrat,GPC.Cr_gamadj,Glc.Cr_gamadj,Glu.Cr_gamadj,GPC.Cho.Cr_gamadj,GABA.Cr_gamadj,NAA.Cr_gamadj,mI.Cr_gamadj,Gln.Cr_gamadj,NAAG.Cr_gamadj,Glu.Gln.Cr_gamadj,Glu.Cr)
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
  luna_data <- luna_data[, c("met", "ld8", "region", "age", "GMrat","Cr")]
  #pivot wide for easy comparison with other spreadsheet
  luna_data <- luna_data %>%
    pivot_wider(names_from = met, values_from = Cr)
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
  
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGlu0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), Glu ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThGABA0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), GABA ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), mI ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThmI0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), mI ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  
  ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGln0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Thalamus'), Gln ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  
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
  Mth <- rbind(ThGlu,ThGABA,ThmI,ThGln,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  #############################
  ####### Caudate #############
  #############################
  
  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGlu0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), Glu ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CaGABA0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), GABA ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), mI ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CamI0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), mI ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  
  CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Gln ~ condition*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGln0 <- tidy(lmerTest::lmer(data = df %>% filter(group != 'HC' & roi == 'Caudate'), Gln ~ group*timepoint + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  
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
  
  Mca <- rbind(CaGlu,CaGABA,CamI,CaGln,CaGluGln,CaNAAG,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  
  # 9 metabolites x 2 regions = N/18 for FDR
  Mall_NR_BL <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_BL')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_R_BL <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_BL')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_NR_FU <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_FU')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_R_FU <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_FU')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_NR_BL_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_BL:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_R_BL_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_BL:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_NR_FU_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionNR_FU:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_R_FU_r <- rbind(Mth,Mca) %>% filter(term %in% c('conditionR_FU:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  
  
  M <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu ~ condition*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(M, ~ condition | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
}
############################
#### High/Low DOI ##########
############################

if (hilowdoi==T){
  
  df <- df %>% mutate(group = case_when(doi_yr >= median(doi_yr,na.rm=TRUE) ~ 'high',
                                          doi_yr < median(doi_yr,na.rm=TRUE) ~ 'low'))
  df$group[is.na(df$group) & df$condition == 'HC_BL'] = 'HC'
  df$group <- relevel(as.factor(df$group),ref='HC')
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), mI ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus'), GPC ~ group*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGln,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  
  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), mI ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu.Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate'), GPC ~ group*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GPC.Cho ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGln,CaGluGln,CaNAAG,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  Mall_high_r <- rbind(Mth1,Mca1) %>% filter(term %in% c('grouphigh:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_high <- rbind(Mth1,Mca1) %>% filter(term %in% c('grouphigh')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_low_r <- rbind(Mth1,Mca1) %>% filter(term %in% c('grouplow:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_low <- rbind(Mth1,Mca1) %>% filter(term %in% c('grouplow')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  
}


#############################
### Group Level ###
#############################

if (group==T){
  
  
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GABA ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), mI ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Gln ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAAG ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), NAA ~ group_level*hemi + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus' & timepoint == 'BL'), GPC ~ group_level*hemi + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group_level + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGln,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  
  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), mI ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  #CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), NAA ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC ~ group_level*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & timepoint == 'BL'), GPC.Cho ~ group_level*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGln,CaGluGln,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  Mall_R_SZ <- rbind(Mth1,Mca1) %>% filter(term %in% c('group_levelSZ:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr) %>% filter(pfdr < 0.05)
  Mall_SZ <- rbind(Mth1,Mca1) %>% filter(term %in% c('group_levelSZ')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value))) %>% arrange(roi,pfdr) %>% filter(pfdr < 0.05)
  
  
  library(emmeans)
  
  # examine emmeans
  MCaGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGlu, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")

  # examine emmeans
  MCaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
    
  # examine emmeans
  MThGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGlu, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGluGln, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
 
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGluGln, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThGluGln, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans, yes there is more in L in HC and more in R in SZ
  # MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  # emm <- emmeans(MThNAAG, ~ hemi | group_level, data = df %>% filter(roi == 'Thalamus'))
  # # Convert to data frame
  # emm_df <- as.data.frame(emm)
  # pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MThNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAA, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
   
  # examine emmeans
  MThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MThNAAG, ~ group_level | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  # examine emmeans
  MCamI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), mI ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCamI, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate'))
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
  MCaGluGln <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu.Gln ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGluGln, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  MCaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ group_level*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(MCaGABA, ~ group_level | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
}

#####################
###.  E/I Balance ###
#####################

if (eibalance){
  
  df <- df %>% group_by(group_level) %>% mutate(GABA_sc = scale(GABA), Glu_sc = scale(Glu))
  
  
  #EI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ scale(Glu.Gln)*group_level*hemi + sex + scale(GMrat) + (1|id))
  EI1 <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate' & group != 'HC'), GABA_sc ~ Glu_sc:group:timepoint + group*timepoint + hemi + age_sc + sex + scale(GMrat) + (1|id))
  #EI2 <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ scale(Gln)*group_level*hemi + sex + scale(GMrat) + (1|id))
  
  
  #EI <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ scale(Glu.Gln)*group + hemi + sex + scale(GMrat) + (1|id))
  EI2 <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus' & group != 'HC'), GABA_sc ~ Glu_sc*group + timepoint + hemi + sex + scale(GMrat) + (1|id))
  #EI2 <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ scale(Gln)*group + hemi + sex + scale(GMrat) + (1|id))
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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
  
  #####################
  #### left side #####
  #####################
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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
  
  
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint','hemi'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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
  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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


  dfw <- df %>% pivot_wider(id_cols = id, names_from = c('group_level','roi','timepoint'), values_from = c('GMrat','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'))
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


  
}

#############################
### Group Level R/NR.     ###
#############################

if (group_r_nr==T){
  
  df$group <- relevel(as.factor(df$group),'HC')
  
  # 2026-06-04 Will need to add hemi eventually
  ThGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), mI ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), Glu.Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  ThGpc <- tidy(lm(data = df %>% filter(roi == 'Thalamus'), GPC ~ group*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  #ThCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GPC.Cho ~ group + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGln,ThGluGln,ThNAAG,ThNAA,ThGpc) %>% mutate(roi = 'Thalamus')
  
  
  CaGlu <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), mI ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu.Gln ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaGpc <- tidy(lm(data = df %>% filter(roi == 'Caudate'), GPC ~ group*hemi + sex + scale(GMrat))) %>% mutate(metabolite = 'GPC')
  CaCho <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GPC.Cho ~ group*hemi + sex + scale(GMrat) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGln,CaGluGln,CaNAAG,CaNAA,CaGpc,CaCho) %>% mutate(roi = 'Caudate')
  
  
  Mall_R_Rem <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupR:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_Rem <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_R_NR <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupNR:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_NR <- rbind(Mth1,Mca1) %>% filter(term %in% c('groupNR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  
  CaGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaGABA, ~ group | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  ThGABA <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), GABA ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(ThGABA, ~ group | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  CaGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaGlu, ~ group | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
 
  CaGlu <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), Glu ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaGlu, ~ hemi | group, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
   
  CamI <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), mI ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CamI, ~ group | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  CaNAA <- lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), NAA ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(CaNAA, ~ group | hemi, data = df %>% filter(roi == 'Caudate'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
  ThNAAG <- lmerTest::lmer(data = df %>% filter(roi == 'Thalamus'), NAAG ~ group*hemi + sex + scale(GMrat) + (1|id))
  emm <- emmeans(ThNAAG, ~ group | hemi, data = df %>% filter(roi == 'Thalamus'))
  # Convert to data frame
  emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  
}

if (clinical==T){
  
  rm(df)
  
  setwd('/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper')
  
  #read in ALL corrected data (including duplicates)
  #isolate patients with baseline and follow-up only (ie only patients with known remitter status)
  sarpal_data_all <- read_csv('sarpal_mrsi_adj_all_20250310.csv') %>%
    filter(source == 'BL.FU')
  
  #add in sex from outside spreadsheet
  # 2026-06-02 AndyP changed to sheet 3 was reading incorrect sheet
  sarpal_data_supplemental <- read_excel('7T.DARES.baselines.xlsx', sheet = 3) %>% group_by(RECID) %>% slice(1) %>% ungroup()
  
  df <- inner_join(sarpal_data_all,sarpal_data_supplemental,by='RECID')
  
  df <- df %>% mutate(sex = as.factor(case_when(sex == 1 ~ 'M',
                                                sex == 2 ~ 'F')))
  
  
  
  ############################
  #add in remitter status from outside spreadsheet
  supplemental_data2 <- read_excel('BPRS_items_MIKE.xlsx')
  supplemental_data2 <- filter(supplemental_data2) %>% mutate(timepoint = case_when(Timepoint == 0 ~ 1,
                                                                                    Timepoint == 6 ~ 2))
  supplemental_data2 <- supplemental_data2 %>%
    mutate(Remitter_Status = ifelse(timepoint == 2 & (CONCDIS4 <= 3 & HALL12 <= 3 & UTC15 <= 3), "R", "NR"))
  
  df <- inner_join(df,supplemental_data2,by=c('RECID','timepoint'))
  
  sd3 <- read_excel('CGI values.xlsx') %>% mutate(timepoint = case_when(TP == 0 ~ 1,
                                                                        TP == 1 ~ 2))
  sd3$timepoint <- as.double(sd3$timepoint)
  
  df <- inner_join(df, sd3, by=c('RECID','timepoint'))
  df <- df %>% mutate(timepoint = case_when(timepoint == 1 ~ 'BL',
                                            timepoint == 2 ~ 'FU'))
  df$timepoint <- relevel(factor(df$timepoint),ref = 'BL')
  df$group <- relevel(factor(df$Remitter_Status),ref = 'R')
  df$sex <- relevel(factor(df$sex),ref = "F")
  df$doi_yr <- df$`DOI (years)`
  df <- df %>% dplyr::select(RECID,hemi,timepoint,age,doi_yr,group,GMrat,GPC.adj,sex,roi,Glc.adj,Glu.adj,GPC.Cho.adj,GABA.adj,NAA.adj,mI.adj,Gln.adj,NAAG.adj,POSSX,Glu.Gln.adj,Severity,Improvement,`Working Memory`,`Verbal Learning`,`Strauss-Carpenter Total`,`Processing Speed`,`Reasoning-Problem Solving`,`Social Cognition`,`MCCB Overall`)
  df <- df %>% dplyr::rename(id = 'RECID',GPC = 'GPC.adj',Glu = 'Glu.adj',GPC.Cho = 'GPC.Cho.adj',Glc = 'Glc.adj',GABA = 'GABA.adj',NAA = 'NAA.adj',mI = 'mI.adj',Gln = 'Gln.adj',NAAG = 'NAAG.adj',Glu.Gln = 'Glu.Gln.adj')

  df <- df %>% mutate(possx_sc = (POSSX), sev_sc = (Severity), imp_sc = (Improvement), 
                      wm_sc = scale(`Working Memory`), mccb_sc = scale(`MCCB Overall`), vl_sc = scale(`Verbal Learning`),
                      ps_sc = scale(`Processing Speed`),sct_sc = scale(`Strauss-Carpenter Total`),rpc_sc = scale(`Reasoning-Problem Solving`),
                      sc_sc = scale(`Social Cognition`))
  
  df <- df %>% dplyr::select(!`Working Memory` & !`Verbal Learning` & !`Strauss-Carpenter Total` & !`Processing Speed` & !`Reasoning-Problem Solving` & !`Social Cognition` & !`MCCB Overall`)
  dfrem <- df %>% dplyr::select(id,sex,timepoint,doi_yr,wm_sc,possx_sc,sev_sc,imp_sc,mccb_sc,vl_sc,sc_sc,ps_sc,sct_sc,rpc_sc)
  dfrem <- dfrem %>% group_by(id,sex,timepoint) %>% slice(1) %>% ungroup()
  dfrem <- dfrem %>% pivot_wider(id_cols = c('id','sex','doi_yr'), names_from = 'timepoint', values_from = c('possx_sc','sev_sc','imp_sc'))
  dfrem <- dfrem %>% mutate(dpossx = (possx_sc_FU - possx_sc_BL)/possx_sc_BL,
                            dsev = (sev_sc_FU - sev_sc_BL),
                            dint = (imp_sc_FU - imp_sc_BL))
  df <- df %>% pivot_wider(id_cols = c('id','roi','hemi','sex'), names_from = 'timepoint', values_from = c('GMrat','Glu','Glc','Glu','GPC.Cho','GABA','NAA','mI','Gln','NAAG','Glu.Gln'),values_fill = 0)
  df <- df %>% mutate(dGlu = Glu_FU - Glu_BL,
                      dGlc = Glc_FU - Glc_BL,
                      dGPC.Cho = GPC.Cho_FU - GPC.Cho_BL,
                      dGABA = GABA_FU - GABA_BL,
                      dNAA = NAA_FU - NAA_BL,
                      dmI = mI_FU - mI_BL,
                      dGln = Gln_FU - Gln_BL,
                      dNAAG = NAAG_FU - NAAG_BL,
                      dGlu.Gln = Glu.Gln_FU - Glu.Gln_BL)
  dfrem <- dfrem %>% dplyr::select(!sex)
  df <- inner_join(df,dfrem,by=c('id'))
  
  df$possx_sc <- scale(df$dpossx)
  df$doi_yr_sc <- scale(Winsorize(df$doi_yr, val = quantile(df$doi_yr,probs = c(0.05,0.95),na.rm=TRUE)))
  
  df1 <- df #%>% filter(doi_yr >= median(doi_yr,na.rm=TRUE))
  df1 <- df1 %>% mutate(group = case_when(doi_yr >= median(doi_yr,na.rm=TRUE) ~ 'high',
                                          doi_yr < median(doi_yr,na.rm=TRUE) ~ 'low'))
  df1$group <- relevel(as.factor(df1$group),ref='low')
  
  df1$curr_var <- df1$group
  
  ThGlu <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dGlu ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  ThGABA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dGABA ~ curr_var*hemi  + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  ThmI <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dmI ~ curr_var*hemi  + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  ThGln <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dGln ~ curr_var*hemi  + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  ThGluGln <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dGlu.Gln ~ curr_var*hemi  + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  ThNAAG <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dNAAG ~ curr_var*hemi  + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  ThNAA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dNAA ~ curr_var*hemi  + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  #ThCho <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Thalamus'), dGPC.Cho ~ curr_var*hemi  + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mth1 <- rbind(ThGlu,ThGABA,ThmI,ThGln,ThGluGln,ThNAAG,ThNAA) %>% mutate(roi = 'Thalamus')
  
  CaGlu <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate'), dGlu ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Glu')
  CaGABA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate'), dGABA ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GABA')
  CamI <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate'), dmI ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'mI')
  #CaGln <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), dGln ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'Gln')
  CaGluGln <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate'), dGlu.Gln ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GluGln')
  #CaNAAG <- tidy(lmerTest::lmer(data = df %>% filter(roi == 'Caudate'), dNAAG ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAAG')
  CaNAA <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate'), dNAA ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'NAA')
  # convergence issues in this lmer model, low variance at the subject level so just use lm
  CaCho <- tidy(lmerTest::lmer(data = df1 %>% filter(roi == 'Caudate'), dGPC.Cho ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU) + (1|id))) %>% dplyr::select(!df & !group & !effect) %>% mutate(metabolite = 'GPC.Cho')
  Mca1 <- rbind(CaGlu,CaGABA,CamI,CaGluGln,CaNAA,CaCho) %>% mutate(roi = 'Caudate')


  Mall <- rbind(Mth1,Mca1) %>% filter(term %in% c('curr_varhigh')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))
  Mall_r <- rbind(Mth1,Mca1) %>% filter(term %in% c('curr_varhigh:hemiR')) %>% mutate(pfdr = p.adjust(p.value,method = 'fdr',n=length(p.value)))

  M <- lm(data = df1 %>% filter(roi == 'Caudate'), dGPC.Cho ~ curr_var*hemi + sex + scale(GMrat_BL) + scale(GMrat_FU))
  emm <- emmeans(M, ~ curr_var | hemi,data = df1 %>% filter(roi == 'Caudate'))
  # # Convert to data frame
  # emm_df <- as.data.frame(emm)
  pairs(emm,adjust = "fdr")
  emtrends(M, specs=pairwise ~ hemi:curr_var, var = 'curr_var', at = c(-1.5,1.5),data = df %>% filter(roi=='Caudate'))
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

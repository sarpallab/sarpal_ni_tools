library(tidyr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(gratia)
library(mgcViz)
library(patchwork)
library(gamm4)
library(readxl)
library(broom)
library(lubridate)

#basedir <- '/ix1/ginger/dsarpal/lab/reorg/projects/20260626-MRSI-Complete'
# macbook
#basedir <- '/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper'
# mac mini
basedir <- '/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Papers/Working_Drafts/Mike_MRSI_Paper'

setwd(basedir)

# Load models
load('20260706-gammodels.Rdata')
gam_models <- M

load('20260706-gamadj-HC.Rdata')


# Load SZ MRSI data
szmet_orig <- readxl::read_xlsx('sarpal_mrsi_original_07062026.xlsx', sheet = 1) %>%
  rename(scan_date = Scan_date) %>% mutate(source = 'orig')
#szmet_orig <- szmet_orig %>% group_by(RECID) %>% separate_wider_delim(cols = region, delim = " ",names = c('hemisphere','roi'), cols_remove = TRUE,too_few = "align_start") %>% ungroup()
szmet_orig <- szmet_orig %>% rename(roi = 'region')


# Match ROI names, create date value
szmet_new <- szmet_orig %>% 
  separate(roi, c('hemi','roi'), ' ', convert=TRUE) %>%
  mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
  mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
  mutate(dateNumeric = as.numeric(as.POSIXct(scan_date, format="%Y-%m-%d")),
         hemi = ifelse(hemi=='left', 'L', ifelse(hemi=='right', 'R', NA))) %>%
  mutate(region = paste0(hemi,' ', roi))

# this did not fix the negative mI values problem 2026-07-17 AndyP
#szmet_new$scan_date[szmet_new$scan_date > "2022-01-01 UTC"] <- "2022-01-01 UTC"

# Grab model based on ROI & metabolite
gam.model <- gam_models$`model_R Thalamus_mI.Cr`
#metdata <- gam_models %>% filter(met == 'GABA', biregion == 'Thalamus') %>% pull(metdata)
metdata <- met_out1 %>% filter(label == 'R Thalamus')
metdata <- metdata %>% mutate(dateNumeric1 = ymd(dateNumeric))
meandate <- mean(metdata$dateNumeric1,na.rm=TRUE)
# Extract necessary data from sz data frame
this.met <- szmet_new %>%
  select(RECID, timepoint, age, hemi, roi, dateNumeric, GMrat, mI.Cre = `mI/Cre`) %>%
  filter(roi == 'Thalamus')
this.met$mI.Cre <- as.numeric(this.met$mI.Cre)


this.met <- this.met %>% rename(mI.Cr = mI.Cre) %>% 
  mutate(mI.Cr = replace_na(mI.Cr,mean(mI.Cr,na.rm=T))) %>%
  mutate(dateNumeric = as.integer(format(as_datetime(dateNumeric, tz = "UTC"), "%Y%m%d"))) 

## first, get residual (difference from expectation given date, GMrat, age)
yhat <- unname(predict(gam.model,this.met))
e <- this.met$mI.Cr - yhat

this.met <- this.met %>% mutate(dateNumeric = meandate)
this.met <- this.met %>% mutate(dateNumeric = as.numeric(gsub("-","",dateNumeric)))
## now, predict @ mean date (but with real age & GMrat) and add back in residual
this.met$mI.Cr.adj <- unname(predict(gam.model, this.met)) + e

ggplot(data = this.met, aes(x = mI.Cr, y = mI.Cr.adj, color=as.factor(hemi), shape=as.factor(timepoint))) +
  geom_point() +
  stat_smooth(method='lm', alpha=0.2) +
  geom_abline(slope=1, intercept = 0, linetype=2) +
  theme(legend.position=c(.1, .8)) +
  theme_bw()


# Loop all
names(szmet_new %>% select(-contains('%') & contains('/Cre')))

sz_mets <- c('Cre','GABA/Cre', 'Glu/Cre', 'Gln/Cre','Glc/Cre', 'GSH/Cre','mI/Cre','NAA/Cre', 'Tau/Cre', 'Glu.Gln/Cre', 'NAAG/Cre', 'GPC/Cre', 'GPC.Cho/Cre')
mod_mets <- c('Cr','GABA.Cr',    'Glu.Cr',     'Gln.Cr',    'Glc.Cr',     'GSH.Cr',    'mI.Cr',    'NAA.Cr',     'Tau.Cr',     'Glu.Gln.Cr', 'NAAG.Cr', 'GPC.Cr','GPC.Cho.Cr')
roiset <- c('L Caudate','R Caudate','L Thalamus','R Thalamus')

#Glu+Gln, NAAG, GPC, GPC+CHO

sz_met_out <- NULL
for (thisroi in roiset) {
  print(thisroi)
  for (meti in seq(1, length(sz_mets))) {
    this_mod_met <- mod_mets[meti]
    print(this_mod_met)
    this_sz_met <- sz_mets[meti]
    
    # Grab model based on ROI & metabolite
    modelstr <- paste0('model_',thisroi,'_',mod_mets[meti],'')
    gam.model <- gam_models[[modelstr]]
    metdata <- met_out1 %>% filter(label == thisroi)
    metdata <- metdata %>% mutate(dateNumeric1 = ymd(dateNumeric))
    meandate <- mean(metdata$dateNumeric1,na.rm=TRUE)
    # Extract necessary data from sz data frame
    this.met <- szmet_new %>%
      select('RECID', 'timepoint', 'age', 'region', 'dateNumeric', 'GMrat', all_of(this_sz_met)) %>%
      filter(region == thisroi)
    this.met[[this_sz_met]] <- as.numeric(this.met[[this_sz_met]])
    
    new_col_name <- this_mod_met
    
    this.met <- this.met %>% 
      mutate(!!new_col_name := replace_na(.data[[this_sz_met]],mean(.data[[this_sz_met]],na.rm=TRUE))) %>%
      mutate(dateNumeric = as.integer(format(as_datetime(dateNumeric, tz = "UTC"), "%Y%m%d"))) 
    
    ## first, get residual (difference from expectation given date, GMrat, age)
    
    tryCatch({
      yhat <- unname(predict(gam.model,this.met))
      e <- this.met[[this_mod_met]] - yhat
      
      this.met <- this.met %>% mutate(dateNumeric = meandate)
      this.met <- this.met %>% mutate(dateNumeric = as.numeric(gsub("-","",dateNumeric)))
      ## now, predict @ mean date (but with real age & GMrat) and add back in residual
      temp <- unname(predict(gam.model, this.met)) + e
      this.met <- this.met %>% mutate(value = temp, metabolite = paste0(this_mod_met,'_gamadj')) %>% select(!all_of(this_sz_met)) %>% select(!all_of(this_mod_met))
      if (length(sz_met_out)==0){
        sz_met_out <- this.met
      } else {
        sz_met_out <- rbind(sz_met_out,this.met)
      }
      
    }, error = function(e) {
      return(NULL)
    })
  }
}

save(sz_met_out, file='20260717-SSD-gamadj.Rdata')

# adj.df.wide <- merge(
#   adj.df %>% select(-met.adj) %>% pivot_wider(names_from = metname, values_from = met),
#   adj.df %>% select(-met) %>% mutate(metname = paste0(metname, '.adj')) %>% pivot_wider(names_from = metname, values_from = met.adj),
#   by = c('RECID','timepoint','age','region','hemi','roi','dateNumeric','GMrat','source'),
#   all=T) %>% 
#   select(RECID, timepoint, age, hemi, roi, dateNumeric, source, GMrat, GABA, GABA.adj, Glu, Glu.adj, Gln, Gln.adj, Glc, Glc.adj, GSH, GSH.adj, mI, mI.adj, NAA, NAA.adj, Tau, Tau.adj)
# 
# dim(adj.df.wide)
# View(adj.df.wide)
# 
# write.csv(x = adj.df.wide, file = '~/scratch/deepak/sarpal_mrsi_adj_all_20250212.csv')
# 
# 
# ggplot(data = adj.df.wide %>% filter(roi == 'Thalamus'), 
#        aes(x = Glu, y = Glu.adj, color=as.factor(hemi), shape=as.factor(timepoint))) +
#   geom_point() +
#   geom_path(aes(group = interaction(RECID, hemi))) +
#   stat_smooth(method='lm', alpha=0.2) +
#   geom_abline(slope=1, intercept = 0, linetype=2) +
#   theme(legend.position=c(.1, .8))
# 
# 
# #gam.model <- gam_models %>% filter(met == 'Glu', biregion == 'DLPFC') %>% pull(model)
# #print(plot(getViz(gam.model[[1]]), allTerms = T), pages = 1)
# 
# head(szmet)
# szmet %>% select(region) %>% distinct()
# 
# unique(gam_models$biregion)
# unique(gam_models$met)
# 
# print(plot(getViz(gam.model[[1]]), allTerms = T), pages = 1)
# 
# 
# szdata <- adj.df.wide %>% filter(roi == 'Caudate')
# gam.model <- gam_models %>% filter(met == 'GABA', biregion == 'Caudate') %>% pull(model)
# metdata <- gam_models %>% filter(met == 'GABA', biregion == 'Caudate') %>% pull(metdata)
# 
# ## adjust original data
# luna_yhat <- unname(predict(gam.model[[1]], metdata[[1]]))
# luna_e <- metdata[[1]]$Cr - luna_yhat
# metdata[[1]]$met.adj <- unname(predict(gam.model[[1]], metdata[[1]] %>% mutate(dateNumeric = mean(metdata[[1]]$dateNumeric, na.rm=T)))) + luna_e
# 
# 
# ## combine
# all.data <- rbind(szdata %>% filter(timepoint == 1) %>% select(age, GABA.adj) %>% mutate(group = 'SZ1'),
#       szdata %>% filter(timepoint == 2) %>% select(age, GABA.adj) %>% mutate(group = 'SZ2'),
#       metdata[[1]] %>% filter(age >= 18) %>% select(age, GABA.adj = met.adj) %>% mutate(group = 'Cn'))
# 
# ## plot
# ggplot(data = all.data, aes(x = age, y = GABA.adj, color = group)) +
#   geom_point() +
#   stat_smooth(method = 'gam', alpha = 0.1)
# 
# ggplot(data = all.data, aes(x = group, y = GABA.adj, fill=group)) + 
# #  geom_violin(draw_quantiles = c(.25, .5, .75))
#   geom_boxplot()
# 
# # Compute the mean and SEM for each group
# data_summary <- all.data %>%
#   group_by(group) %>%
#   summarise(
#     mean_GABA = mean(GABA.adj, na.rm = TRUE),
#     sem_GABA = sd(GABA.adj, na.rm = TRUE) / sqrt(n())
#   )
# 
# # Create the bar plot with error bars
# ggplot(data_summary, aes(x = group, y = mean_GABA, fill = group)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_errorbar(aes(ymin = mean_GABA - sem_GABA, ymax = mean_GABA + sem_GABA), 
#                 width = 0.2) +
#   scale_fill_manual(values = c("SZ1" = "blue", "SZ2" = "red", "Cn" = "gray")) +
#   labs(x = "Group", y = "Glu (Mean ± SEM)") +
#   theme_minimal() +
#   theme(legend.position = "none")
# 

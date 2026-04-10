# 2026-02-24 AndyP
# Script to look at Samira's DWI

library(readxl)
library(broom)
library(tidyverse)
library(DescTools)

df <- read_xlsx('/Users/andrew/Downloads/diff_MK_BG.xlsx')

df$age_sc <- scale(df$Age)
df$perc_change_total_sc <- scale(df$perc_change_total)
df$CLZ_dose_FU <- scale(df$CLZ_dose_FU)
df <- df %>% mutate(Other_AP_Factor = case_when(Other_AP == 0 ~ 'None',
                                                Other_AP > 0 ~ 'Other_AP_Present'))
df$Other_AP_Factor <- relevel(as.factor(df$Other_AP_Factor),ref = 'None')

df$vCa_L_percent <- Winsorize(df$vCa_L_percent)
df$vCa_R <- Winsorize(df$vCa_R)
df$GP_L_percent <- Winsorize(df$GP_L_percent)
df$GP_R <- Winsorize(df$GP_R)
df$NAC_L <- Winsorize(df$NAC_L)
df$NAC_R <- Winsorize(df$NAC_R)
df$vmPu_L <- Winsorize(df$vmPu_L)
df$vmPu_R <- Winsorize(df$vmPu_R)
df$dCa_L <- Winsorize(df$dCa_L)
df$dCa_R <- Winsorize(df$dCa_R)
df$dlpu_L_percent <- Winsorize(df$dlpu_L_percent)
df$dlPu_R <-  Winsorize(df$dlPu_R)

m1 <- tidy(lm(data = df, vCa_L_percent ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'vCa_L_percent') %>% filter(term == 'perc_change_total_sc')
m2 <- tidy(lm(data = df, vCa_R ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'vCa_R') %>% filter(term == 'perc_change_total_sc')
m3 <- tidy(lm(data = df, GP_L_percent ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'GP_L_percent') %>% filter(term == 'perc_change_total_sc')
m4 <- tidy(lm(data = df, GP_R ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'GP_R') %>% filter(term == 'perc_change_total_sc')
m5 <- tidy(lm(data = df, NAC_L ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'NAC_L') %>% filter(term == 'perc_change_total_sc')
m6 <- tidy(lm(data = df, NAC_R ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'NAC_R') %>% filter(term == 'perc_change_total_sc')
m7 <- tidy(lm(data = df, vmPu_L ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'vmPu_L') %>% filter(term == 'perc_change_total_sc')
m8 <- tidy(lm(data = df, vmPu_R ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'vmPu_R') %>% filter(term == 'perc_change_total_sc')
m9 <- tidy(lm(data = df, dCa_L ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'dCa_L') %>% filter(term == 'perc_change_total_sc')
m10 <- tidy(lm(data = df, dCa_R ~ perc_change_total_sc  + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'dCa_R') %>% filter(term == 'perc_change_total_sc')
m11 <- tidy(lm(data = df, dlpu_L_percent ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor )) %>%
  mutate(outcome = 'dlpu_L_percent') %>% filter(term == 'perc_change_total_sc')
m12 <- tidy(lm(data = df, dlPu_R ~ perc_change_total_sc + age_sc + Sex + CLZ_dose_FU +Other_AP_Factor  )) %>%
  mutate(outcome = 'dlPu_R') %>% filter(term == 'perc_change_total_sc')

mall <- rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)
mall <- mall %>% mutate(pfdr = p.adjust(p.value, 'fdr'))

dflong <- df %>% pivot_longer(cols = 2:13)
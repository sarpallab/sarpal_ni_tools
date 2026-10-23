# 2026-10-23 AndyP
# parse MARS scanning directories and save as interpretable format

library(tidyverse)
library(stringr)

df <- read.csv('/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Data/MARS/directories.csv',header=FALSE)

# Split into 3 distinct columns
df_split <- df %>%
  separate_wider_delim(V1, delim = "/", names = c("dir1", "dir2", "dir3","dir4","dir5","dir6","dir7","dir8","id","id_date","scans"),
                       too_few = "align_start",
                       too_many = "drop")

df_split <- df_split %>% 
  filter(!(id %in% c('Pilot_Prisma2','9279Pilot'))) %>% 
  filter(!is.na(id_date)) %>%
  # 1. Extract the patient ID
  mutate(patient_id = str_remove(id_date, "_.*")) %>%
  
  # 2. Group by each patient
  group_by(patient_id) %>%
  
  # 3. Sort by the id_date string within the group
  arrange(id_date, .by_group = TRUE) %>%
  
  # 4. Assign identical sequential ranks across duplicate rows
  mutate(scan_number = dense_rank(id_date)) %>%
  
  # 5. Clean up grouping
  ungroup()

# Assuming your data frame is named df_split and the target column is 'x'
dfsummary <- df_split %>% 
  mutate(x = scans) %>%
  group_by(patient_id, scan_number) %>%
  summarize(
    # 1. "T1" or "T1w"
    has_T1_or_T1w = any(str_detect(x, "T1|T1w"), na.rm = TRUE),
    
    # 2. "rsfMRI"
    has_rsfMRI = any(str_detect(x, "rsfMRI"), na.rm = TRUE),
    
    # 3. "RISEAssociativeRecognition"
    has_RISEAssociativeRecognition = any(str_detect(x, "RISEAssociativeRecognition"), na.rm = TRUE),
    
    # 4. "sL_dA_SUM"
    has_sL_dA_SUM = any(str_detect(x, "sL_dA_SUM"), na.rm = TRUE),
    
    # 5. "CPT5" or "CPT4"
    has_CPT5_or_CPT4 = any(str_detect(x, "CPT5|CPT4"), na.rm = TRUE),
    
    # 6. "SpinEchoFieldMap_PA"
    has_SpinEchoFieldMap_PA = any(str_detect(x, "SpinEchoFieldMap_PA"), na.rm = TRUE),
    
    # 7. "SpinEchoFieldMap_AP"
    has_SpinEchoFieldMap_AP = any(str_detect(x, "SpinEchoFieldMap_AP"), na.rm = TRUE),
    
    # 8. "RISEAssociativeRecognition" and "PhysioLog"
    has_RISEAssociative_and_Physio = any(str_detect(x, "RISEAssociativeRecognition") & str_detect(x, "PhysioLog"), na.rm = TRUE),
    
    # 9. "AX-CPT"
    has_AX_CPT = any(str_detect(x, "AX-CPT"), na.rm = TRUE),
    
    # 10. "rsfMRI" and "PhysioLog"
    has_rsfMRI_and_Physio = any(str_detect(x, "rsfMRI") & str_detect(x, "PhysioLog"), na.rm = TRUE),
    
    # 11. "AX-CPT" and "PhysioLog"
    has_AX_CPT_and_Physio = any(str_detect(x, "AX-CPT") & str_detect(x, "PhysioLog"), na.rm = TRUE),
    
    # 12. ("CPT5" or "CPT4") and "PhysioLog"
    has_CPT_and_Physio = any(str_detect(x, "CPT5|CPT4") & str_detect(x, "PhysioLog"), na.rm = TRUE),
    
    # 13. "x_sL_dA"
    has_x_sL_dA = any(str_detect(x, "x_sL_dA"), na.rm = TRUE),
    
    # 14. "x_sL_hp"
    has_x_sL_hp = any(str_detect(x, "x_sL_hp"), na.rm = TRUE),
    
    # 15. "RISEImage"
    has_RISEImage = any(str_detect(x, "RISEImage"), na.rm = TRUE),
    
    # 16. "RISEImage" and "PhysioLog"
    has_RISEImage_and_Physio = any(str_detect(x, "RISEImage") & str_detect(x, "PhysioLog"), na.rm = TRUE),
    
    # 17. "RISEItemRecognition"
    has_RISEItemRecognition = any(str_detect(x, "RISEItemRecognition"), na.rm = TRUE),
    
    # 18. "RISEItemRecognition" and "PhysioLog"
    has_RISEItemRecognition_and_Physio = any(str_detect(x, "RISEItemRecognition") & str_detect(x, "PhysioLog"), na.rm = TRUE),
    
    .groups = "drop" # Natively ungroups your data frame so you don't need a separate %>% ungroup()
  )


setwd('/Users/andrew/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SARPALlab - Documents/Data/MARS')
write.csv(dfsummary, file = '2026-10-23-MARS-scan-inventory.csv')
write.csv(df_split, file = '2026-10-23-MARS-scan-raw-folder-list.csv')

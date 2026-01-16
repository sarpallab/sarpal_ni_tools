# 2026-01-16 AndyP
# Generate GUID mapping to subject for bids2nda
# format apparently needs to be subject - NDARXXXX, see https://github.com/WillForan/bids2nda/blob/lncd/examples/guid_map.txt

library(tidyverse)

ids <- read.csv('/Users/andrew/Desktop/2026_Enact_GUIDs.csv',header=TRUE)

ids$idsout <- apply(ids[, c("src_subject_id", "subjectkey")], 1, function(row) {
  paste(row["src_subject_id"], row["subjectkey"], sep = " - ")
})

idsout <- data.frame(ids = ids$idsout)
write.table(idsout, file='2026-ENACT-GUIDs.txt',sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

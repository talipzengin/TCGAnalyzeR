# Save rda file for Gene ID conversion using BioMART
wd <- paste0("DEA_ref")
dir.create(wd)
setwd(wd)

pd <- dirname(getwd())

## ID conversion

library("biomaRt")
httr::set_config(httr::config(ssl_verifypeer = FALSE))

human_ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
save(human_ensembl, file = "human_ensembl.rda")

print(paste0("########### ", "DEA_ref HAS BEEN PREPARED !!!", " ###########"))

save.image()
q(save="yes")

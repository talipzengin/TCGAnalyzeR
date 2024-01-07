library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least one arguments (TCGA cancer type) must be supplied.n", call.=FALSE)
}

cancer <- args[1]
project <- paste0("TCGA-", cancer)
project_ <- gsub("-","_", project)

data_dir <- "/media/data02/GDC_data/TCGA/"

wd <- paste0(project, "/clinical")
dir.create(wd)
setwd(wd)

load(file = paste0(data_dir, project, "/", project_, "_clinical_data.rda"))
  
clinic <- clinic[order(clinic$submitter_id),]

clinic[clinic == "Not Reported"] <- NA
clinic[clinic == "not reported"] <- NA

all_na <- function(x) any(!is.na(x))
cfu <- clinic %>% select_if(all_na)

cfu <- cfu %>% mutate_if(is.character, as.factor)
save(cfu, file = paste0(cancer, "_parsed_clinical_data.rda"))
write.table(cfu, file = paste0(cancer, "_parsed_clinical_data.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

print(paste0("####### ", cancer, " CLINICAL ANALYSIS COMPLETED", " #######"))

save.image()
q(save="yes")

library(TCGAbiolinks)
library(gaia)
library(downloader)
library(readr)
library(data.table)
library(statmod)
library(tidyverse)

wd <- paste0("CNV_ref")
dir.create(wd)
setwd(wd)

pd <- dirname(getwd())

########################################################
download.file("https://api.gdc.cancer.gov/data/77fbfff6-2acc-47ca-a5f6-c488beb46879", destfile="snp6.na35.liftoverhg38.txt.zip")
file <- unzip("snp6.na35.liftoverhg38.txt.zip")
markersMatrix <- readr::read_tsv(file, col_names = TRUE, col_types = "ccn____", progress = TRUE)
colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")
unique(markersMatrix$Chromosome)
markersMatrix[markersMatrix$Chromosome == "X","Chromosome"] <- "23"
markersMatrix[markersMatrix$Chromosome == "Y","Chromosome"] <- "24"
markersMatrix$Chromosome <- as.integer(markersMatrix$Chromosome)
markersMatrix <- markersMatrix[order(markersMatrix$Chromosome, markersMatrix$Start),]
markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")
markersMatrix <- markersMatrix[!duplicated(markerID),] # Removed duplicates
rm(markerID)
save(markersMatrix, file = paste0(pd, "/CNV_ref/markersMatrix.rda"), compress = "xz")
write.table(markersMatrix, file = paste0(pd, "/CNV_ref/markersMatrix.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
#######################################################

## Gene annotation of recurrent CNV

library(GenomicRanges)

####### Get gene information from GENCODE using biomart #############
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38") 
genes <- genes[genes$external_gene_name != "" & genes$chromosome_name %in% c(1:22,"X","Y"),]
genes[genes$chromosome_name == "X", "chromosome_name"] <- 23
genes[genes$chromosome_name == "Y", "chromosome_name"] <- 24
genes$chromosome_name <- sapply(genes$chromosome_name,as.integer)
genes <- genes[order(genes$start_position),]
genes <- genes[order(genes$chromosome_name),]
genes <- genes[,c("external_gene_name", "chromosome_name", "start_position","end_position", "entrezgene_id")]
colnames(genes) <- c("GeneSymbol","Chr","Start","End", "EntrezGeneID")
genes_GR <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
save(genes_GR, genes, file = paste0(pd, "/CNV_ref/genes_GR.rda"), compress = "xz")
#################################################################################

print(paste0("########### ", "CNV_ref HAS BEEN PREPARED !!!", " ###########"))

q(save="yes")

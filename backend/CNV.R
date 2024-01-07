## CNV Analysis of TCGA ----------------------------------------------------

library(TCGAbiolinks)
library(gaia)
library(downloader)
library(readr)
library(data.table)
library(statmod)
library(tidyverse)

load_cnv <- function (segmentation_matrix, markers_list, num_of_samples) {
  message("Loading Copy Number Data")
  chromosomes <- as.numeric(sort(unique(names(markers_list))))
  chromosomes <- chromosomes[which(!is.na(chromosomes))]
  aberration_kinds <- 1:length(unique(segmentation_matrix[,6]))
  names(aberration_kinds) <- sort(unique(segmentation_matrix[,6]))
  samples <- 1:num_of_samples
  if (!is.numeric(segmentation_matrix[, 1])) {
    sample_names <- unique(segmentation_matrix[, 1])
    for (i in 1:length(sample_names)) {
      segmentation_matrix[which(segmentation_matrix[, 1] == 
                                  sample_names[i]), 1] <- i
    }
  }
  region_list <- list()
  for (k in 1:length(aberration_kinds)) {
    region_list[[aberration_kinds[k]]] <- list()
    for (i in 1:length(chromosomes)) {
      region_list[[aberration_kinds[k]]][[chromosomes[i]]] <- matrix(0, 
                                                                     length(samples), ncol(markers_list[[chromosomes[i]]]))
      rownames(region_list[[aberration_kinds[k]]][[chromosomes[i]]]) <- samples
      colnames(region_list[[aberration_kinds[k]]][[chromosomes[i]]]) <- c(1:ncol(markers_list[[chromosomes[i]]]))
    }
  }
  for (k in 1:length(aberration_kinds)) {
    ab_ids <- which(segmentation_matrix[, 6] == names(aberration_kinds[k]))
    tmp_matrix1 <- segmentation_matrix[ab_ids, ]
    if (class(tmp_matrix1) == "numeric") {
      tmp_matrix1 <- t(as.matrix(tmp_matrix1))
    }
    for (i in 1:length(chromosomes)) {
      tmp_matrix2 <- tmp_matrix1[which(tmp_matrix1[, 2] == 
                                         chromosomes[i]), ]
      if (class(tmp_matrix2) == "numeric") {
        tmp_matrix2 <- t(as.matrix(tmp_matrix2))
      }
      message(".", appendLF = FALSE)
      for (j in 1:length(samples)) {
        tmp_matrix3 <- tmp_matrix2[which(tmp_matrix2[, 
                                                     1] == samples[j]), ]
        if (class(tmp_matrix3) == "numeric") {
          tmp_matrix3 <- t(as.matrix(tmp_matrix3))
        }
        if (nrow(tmp_matrix3) > 0) {
          for (t in 1:nrow(tmp_matrix3)) {
            start_prob <- tmp_matrix3[t, 3]
            end_prob <- tmp_matrix3[t, 4]
            start_index <- which(markers_list[[chromosomes[i]]][1,] == start_prob)
            end_index <- which(markers_list[[chromosomes[i]]][2,] == end_prob)
            if ((length(start_index) > 0 ) & (length(end_index) > 0 )){
              region_list[[aberration_kinds[k]]][[chromosomes[i]]][samples[j], 
                                                                   start_index:end_index] <- 1
              
            }
          }
        }
      }
    }
  }
  message("\nDone")
  names(region_list) <- names(aberration_kinds)
  return(region_list)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least one arguments (TCGA cancer type) must be supplied.n", call.=FALSE)
}

cancer <- args[1]
project <- paste0("TCGA-", cancer)
project_ <- gsub("-","_", project)

data_dir <- "/media/data02/GDC_data/TCGA/"

wd <- paste0(project, "/CNV")
dir.create(wd)
setwd(wd)

pd <- dirname(getwd())
ppd <- dirname(pd)

load(file = paste0(data_dir, project, "/", project_, "_Copy_Number_Variation_Masked_Copy_Number_Segment", ".rda"))

cases_cnv <- unique(data$Sample)

# Filter TCGA Replicate Samples (cases_cnv -> uniq_tsb_cnv)
source(paste0(ppd,"/tcga_replicateFilter.R"))
write.table(uniq_tsb_cnv, file = paste0(cancer, "_cases_cnv.tsv"), row.names=FALSE, na="NA", col.names=FALSE, sep="\t")

data <- data %>% filter(Sample %in% uniq_tsb_cnv)

save.image()

# Prepare CNV matrix
cnvMatrix <- data
patients <- unique(substr(cnvMatrix$Sample, 1, 12))
patients <- sort(patients)

# Add label (0 for loss, 1 for gain)
cnvMatrix <- cbind(cnvMatrix, Label=NA)
cnvMatrix[cnvMatrix[,"Segment_Mean"] < -0.3,"Label"] <- 0
cnvMatrix[cnvMatrix[,"Segment_Mean"] > 0.3,"Label"] <- 1
cnvMatrix <- cnvMatrix[!is.na(cnvMatrix$Label),]

# Remove "GDC_Aliquot" and "Segment_Mean" columns and change col.names with reordering
cnvMatrix <- cnvMatrix[,-c(1, 6)]
colnames(cnvMatrix) <- c("Chromosome", "Start", "End", "Num.of.Markers", "Sample.Name", "Aberration")
cnvMatrix <- cnvMatrix[, c(5, 1, 2, 3, 4, 6)]

# Substitute Chromosomes "X" and "Y" with "23" and "24"
cnvMatrix[cnvMatrix$Chromosome == "X","Chromosome"] <- 23
cnvMatrix[cnvMatrix$Chromosome == "Y","Chromosome"] <- 24
cnvMatrix$Chromosome <- as.integer(cnvMatrix$Chromosome)
save(cnvMatrix, file=paste0(cancer, "_cnvMatrix.rda"))
write.table(cnvMatrix, file = paste0(cancer, "_cnvMatrix.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

save.image()

########################################################
#download.file("https://api.gdc.cancer.gov/data/77fbfff6-2acc-47ca-a5f6-c488beb46879", destfile="snp6.na35.liftoverhg38.txt.zip")
#file <- unzip("snp6.na35.liftoverhg38.txt.zip")
#markersMatrix <- readr::read_tsv(file, col_names = TRUE, col_types = "ccn____", progress = TRUE)
#colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")
#unique(markersMatrix$Chromosome)
#markersMatrix[markersMatrix$Chromosome == "X","Chromosome"] <- "23"
#markersMatrix[markersMatrix$Chromosome == "Y","Chromosome"] <- "24"
#markersMatrix$Chromosome <- as.integer(markersMatrix$Chromosome)
#markersMatrix <- markersMatrix[order(markersMatrix$Chromosome, markersMatrix$Start),]
#markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")
#markersMatrix <- markersMatrix[!duplicated(markerID),] # Removed duplicates
#rm(markerID)
#save(markersMatrix, file = paste0(ppd, "/CNV_ref/markersMatrix.rda"), compress = "xz")
#write.table(markersMatrix, file = paste0(ppd, "/CNV_ref/markersMatrix.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
#######################################################

load(paste0(ppd, "/CNV_ref/markersMatrix.rda"))

library(gaia)

set.seed(1)
markers_obj <- load_markers(as.data.frame(markersMatrix))
nbsamples <- length(unique(cnvMatrix$Sample))
cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
suppressWarnings(results <- runGAIA(cnv_obj,
                                    markers_obj,
                                    output_file_name = paste0(cancer, "_gaia.txt"),
                                    aberrations = -1,  # -1 to all aberrations
                                    chromosomes = -1, # -1 to all chromosomes
                                    approximation = TRUE, # Set to TRUE to speed up the time requirements
                                    num_iterations = 5000, # Reduced to speed up the time requirements
                                    threshold = 0.05)
)

# Set q-value threshold
threshold <- 0.05       #arg input from user in Rshiny !!!

# Plot the results
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score = 0)
minval <- format(min(RecCNV[RecCNV[,"q-value"] != 0,"q-value"]), scientific = FALSE)
minval <- substring(minval,1, nchar(minval) - 1)
RecCNV[RecCNV[,"q-value"] == 0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"] == as.numeric(minval),]

# Filter results according to threshold
sCNV <- RecCNV[RecCNV[,"q-value"] <= threshold,c(1:4,6)]
sCNV <- sCNV[order(sCNV[,3]),]
sCNV <- sCNV[order(sCNV[,1]),]
colnames(sCNV) <- c("Chr","Aberration","Start","End","q-value")

png(file = paste0(cancer, "_CNVplot.png"), width = 4, height = 4, units = 'in', res = 300)
gaiaCNVplot(RecCNV,threshold)
dev.off()

save(results, RecCNV, threshold, sCNV, file = paste0(cancer, "_CNV_results.rda"))

save.image()

## Gene annotation of recurrent CNV

library(GenomicRanges)

####### Get gene information from GENCODE using biomart #############
#genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38") 
#genes <- genes[genes$external_gene_name != "" & genes$chromosome_name %in% c(1:22,"X","Y"),]
#genes[genes$chromosome_name == "X", "chromosome_name"] <- 23
#genes[genes$chromosome_name == "Y", "chromosome_name"] <- 24
#genes$chromosome_name <- sapply(genes$chromosome_name,as.integer)
#genes <- genes[order(genes$start_position),]
#genes <- genes[order(genes$chromosome_name),]
#genes <- genes[,c("external_gene_name", "chromosome_name", "start_position","end_position", "entrezgene_id")]
#colnames(genes) <- c("GeneSymbol","Chr","Start","End", "EntrezGeneID")
#genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
#save(genes_GR, genes, file = paste0(ppd, "/CNV_ref/genes_GR.rda"), compress = "xz")
#################################################################################

load(paste0(ppd, "/CNV_ref/genes_GR.rda"))

## Recurrent CNV annotation ## 

# Get gene information from GENCODE using biomart data(genes_GR)

sCNV_GR <- makeGRangesFromDataFrame(sCNV, keep.extra.columns = TRUE)

hits <- findOverlaps(genes_GR, sCNV_GR, type = "within")
sCNV_ann <- cbind(sCNV[subjectHits(hits),],genes[queryHits(hits),])
AberrantRegion <- paste0(sCNV_ann[,1],":",sCNV_ann[,3],"-",sCNV_ann[,4])
GeneRegion <- paste0(sCNV_ann[,7],":",sCNV_ann[,8],"-",sCNV_ann[,9])
AmpDel_genes <- cbind(sCNV_ann[,c(6,2,5)],AberrantRegion,GeneRegion,sCNV_ann[,10])
colnames(AmpDel_genes)[6] <- "EntrezGeneID"
AmpDel_genes[AmpDel_genes[,2] == 0,2] <- "Del"
AmpDel_genes[AmpDel_genes[,2] == 1,2] <- "Amp"
rownames(AmpDel_genes) <- NULL
write.table(AmpDel_genes, file = paste0(cancer, "_AmpDel_genes.tsv"), sep = "\t", col.names = TRUE, row.names = TRUE)

Amp_genes <- AmpDel_genes[AmpDel_genes[,2] == "Amp",]
write.table(Amp_genes, file = paste0(cancer, "_Amp_genes.tsv"), sep = "\t", col.names = TRUE, row.names = TRUE)
Del_genes <- AmpDel_genes[AmpDel_genes[,2] == "Del",]
write.table(Del_genes, file = paste0(cancer, "_Del_genes.tsv"), sep = "\t", col.names = TRUE, row.names = TRUE)

save(AmpDel_genes, Amp_genes, Del_genes, file = paste0(cancer, "_CNV_results_genes.rda"))

save.image()

### Enrichment Analysis

library(clusterProfiler)
library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)

AmpDel_genes2 <- na.omit(AmpDel_genes)
Amp_genes2 <- AmpDel_genes2[AmpDel_genes2[,2] == "Amp",]
Del_genes2 <- AmpDel_genes2[AmpDel_genes2[,2] == "Del",]

## KEGG over-representation test
EA <- enrichKEGG(gene = AmpDel_genes2$EntrezGeneID,
                 keyType = 'ncbi-geneid',
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
# Input ID type can be 'kegg', 'ncbi-geneid', 'ncbi-proteinid' or 'uniprot'.

EA <- setReadable(EA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA), file = paste0(cancer, "_pathways_total.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")


EA_Amp <- enrichKEGG(gene = Amp_genes2$EntrezGeneID,
                    keyType = 'ncbi-geneid',
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

EA_Amp <- setReadable(EA_Amp, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA_Amp), file = paste0(cancer, "_pathways_Amp.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")



EA_Del <- enrichKEGG(gene = Del_genes2$EntrezGeneID,
                      keyType = 'ncbi-geneid',
                      organism = 'hsa',
                      pvalueCutoff = 0.05)

EA_Del <- setReadable(EA_Del, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA_Del), file = paste0(cancer, "_pathways_Del.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

save.image()

## Patient based CNVs of significant genes
cnvMatrix_ <- cnvMatrix[order(cnvMatrix$Chromosome, cnvMatrix$Start),]

cnv_results <- AmpDel_genes2
cnv_results <- cnv_results[!duplicated(cnv_results$GeneSymbol),]
cnv_results <- cnv_results[order(cnv_results$`q-value`),]

nrow <- length(patients)
ncol <- length(unique(cnv_results$GeneSymbol))
cnv_gene_matrix <- data.frame(matrix(data="", ncol = ncol, nrow = nrow), stringsAsFactors=FALSE)
colnames(cnv_gene_matrix) <- unique(cnv_results$GeneSymbol)
rownames(cnv_gene_matrix) <- patients

for (i in 1:nrow(cnv_results)){
  string <- cnv_results$GeneRegion[i]
  pattern <- "[[:punct:]]"
  splits <- strsplit(string, pattern)
  Chr <- splits[[1]][1]
  Chr_start <- splits[[1]][2]
  Chr_end <- splits[[1]][3]
  gene_name <- cnv_results$GeneSymbol[i]
  for (j in 1:nrow(cnvMatrix_)){
    if (as.numeric(Chr) == cnvMatrix_$Chromosome[j] && as.numeric(Chr_start) >= cnvMatrix_$Start[j] && as.numeric(Chr_end) <= cnvMatrix_$End[j]){
      patient <- substr(cnvMatrix_$Sample.Name[j], 1, 12)
      if(cnvMatrix_$Aberration[j] == 0){
        gene_abr <- "Del;"
      } else {
        gene_abr <- "Amp;"
      }
      cnv_gene_matrix[patient,gene_name] <- gene_abr
    }
  }
}

save(cnv_gene_matrix, file = paste0(cancer, "_cnv_gene_matrix.rda"))
write.table(cnv_gene_matrix, file = paste0(cancer, "_cnv_gene_matrix.tsv"), row.names=TRUE, na="NA", col.names=TRUE, sep="\t")

save.image()

print(paste0("########## ", cancer, " CNV ANALYSIS COMPLETED", " ##########"))

save.image()
q(save="yes")

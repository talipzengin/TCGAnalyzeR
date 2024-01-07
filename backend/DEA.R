## GENE EXPRESSION DATA ANALYSIS --------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two arguments (TCGA cancer type & cohort type) must be supplied.n", call.=FALSE)
}

cancer <- args[1]
cohort <- args[2]
project <- paste0("TCGA-", cancer)
project_ <- gsub("-","_", project)

data_dir <- "/media/data02/GDC_data/TCGA/"

wd <- paste0(project, "/DEA_", cohort)
dir.create(wd)
setwd(wd)

pd <- dirname(getwd())
ppd <- dirname(pd)

load(file = paste0(data_dir, project, "/", project_, "_Gene_Expression_Quantification_HTSeq_Counts", ".rda"))

cases_exp <- colData(data)$barcode

if (cancer == "LAML"){
  tumor_sample_type="TB"
  normal_sample_type="NB"
} else {
  tumor_sample_type="TP"
  normal_sample_type="NT"
}

if (cohort == "total"){
  cases_exp <- TCGAquery_SampleTypes(cases_exp, c(normal_sample_type, tumor_sample_type))
} else {
  cases_exp <- TCGAquery_MatchedCoupledSampleTypes(cases_exp, c(normal_sample_type, tumor_sample_type))
}

# Filter TCGA Replicate Samples (cases_exp -> uniq_tsb_exp)
source(paste0(ppd,"/tcga_replicateFilter.R"))
write.table(uniq_tsb_exp, file = paste0(cancer, "_cases_exp.tsv"), row.names=FALSE, na="NA", col.names=FALSE, sep="\t")


data <- subset(data, select = colData(data)$barcode %in% uniq_tsb_exp)


GeneExp_Matrix <- assay(data)
save(GeneExp_Matrix, file = paste0(cancer, "_GeneExp_Matrix.rda"))
write.table(GeneExp_Matrix, file= paste0(cancer, "_GeneExp_Matrix.tsv"), 
            row.names=TRUE, na="NA",col.names=TRUE, sep="\t")

colData <- as.data.frame(colData(data))
write.table(colData, file = paste0(cancer, "_GeneExp_colData.tsv"),
            row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

write.table(colData[,2:4], file = paste0(cancer, "_targets.tsv"),
            row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

save.image()

### Normalization of all Gene Expression Matrix

library(limma)
library(edgeR)

gem <- DGEList(counts=GeneExp_Matrix)

#TMM normalization
gem_norm <- calcNormFactors(gem, method="TMM")
save(gem_norm, file = paste0(cancer, "_gem_normalized.rda"))
write.table(gem_norm$counts, file = paste0(cancer, "_gem_normalized.tsv"), 
            row.names=TRUE, na="",col.names=TRUE, sep="\t")

gem_norm_log <- log2(gem_norm$counts+1.0001)
write.table(gem_norm_log, file= paste0(cancer, "_", cohort, "_gem_norm_log2.tsv"), row.names=TRUE, na="NA",col.names=TRUE, sep="\t")

save.image()

### Differentially Expression Analysis (DEA) by limma-voom

library(limma)
library(edgeR)

targets <- read.delim(paste0(cancer, "_targets.tsv"), sep = "\t")

dge <- DGEList(counts=GeneExp_Matrix)

#filtering to remove rows that consistently have zero or very low counts
keep <- filterByExpr(dge)
dge <- dge[keep,]

#TMM normalization
dgen <- calcNormFactors(dge, method="TMM")
save(dgen, file = paste0(cancer, "_dge_normalized.rda"))
write.table(dgen$counts, file = paste0(cancer, "_dge_normalized.tsv"), 
            row.names=TRUE, na="",col.names=TRUE, sep="\t")

save.image()

# Differential Expression Analysis
design <- model.matrix(~ targets$shortLetterCode)

v <- voom(dgen, design, plot=FALSE)

cor <- duplicateCorrelation(v, design, block = targets$patient)
cor$consensus

v <- voom(dgen, design, block = targets$patient, correlation = cor$consensus)

fit <- lmFit(v, design, block = targets$patient, correlation = cor$consensus)
fit <- eBayes(fit)

degs <- topTable(fit, coef=ncol(design), number=Inf, sort.by="P")
degs <- tibble::rownames_to_column(degs, "ensembl_gene_id")
write.table(degs, file = paste0(cancer, "_degs.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

sign_degs <- topTable(fit, coef=ncol(design), number=Inf, sort.by="P", p.value=0.01, lfc=1)
sign_degs <- tibble::rownames_to_column(sign_degs, "ensembl_gene_id")
save(sign_degs, file = paste0(cancer, "_sign_degs.rda"))
write.table(sign_degs, file = paste0(cancer, "_sign_degs.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

save.image()


## Volcano Plot

TCGAVisualize_volcano(degs$logFC, degs$adj.P.Val,
                      filename = paste0(cancer, "_VolcanoPlot.png"),
                      x.cut = 1,
                      y.cut = 0.01,
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " RNA expression fold change (Log2)",
                      legend = "State",
                      title = "Differential Gene Expression",
                      width = 10,
                      height = 6,
                      dpi = 300)


## ID conversion

library("biomaRt")
#httr::set_config(httr::config(ssl_verifypeer = FALSE))

#human_ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
load(paste0(ppd, "/DEA_ref/human_ensembl.rda"))

#all_genes_ids <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'description'), 
#                      filters = 'ensembl_gene_id', 
#                      values = rownames(gem_norm[["counts"]]), 
#                      mart = human_ensembl)

#write.table(all_genes_ids, file = paste0(cancer, "_all_genes_ids.tsv"), 
#            row.names=FALSE, col.names=TRUE, na="NA", sep="\t")

if (length(rownames(sign_degs)) != 0){
sign_degs_ids <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'description'), 
                      filters = 'ensembl_gene_id', 
                      values = sign_degs$ensembl_gene_id, 
                      mart = human_ensembl)

sign_degs_ids <- sign_degs_ids[which(!duplicated(sign_degs_ids$ensembl_gene_id)),]
write.table(sign_degs_ids, file = paste0(cancer, "_sign_degs_ids.tsv"), 
            row.names=FALSE, col.names=TRUE, na="NA", sep="\t")

sign_degs_with_ids <- inner_join(sign_degs, sign_degs_ids, by = "ensembl_gene_id") 
write.table(sign_degs_with_ids, file = paste0(cancer, "_sign_degs_with_ids.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")
}

save.image()

####

if (exists("sign_degs_ids")){
# Expression Matrix
exp_matrix <- gem_norm[["counts"]][sign_degs_ids$ensembl_gene_id,]

# Exp matrix with Gene Symbols
for (i in 1:nrow(exp_matrix)) {
  if (!identical(sign_degs_ids$external_gene_name[which(sign_degs_ids$ensembl_gene_id == rownames(exp_matrix)[i])], "")) {
  rownames(exp_matrix)[i] <- sign_degs_ids$external_gene_name[which(sign_degs_ids$ensembl_gene_id == rownames(exp_matrix)[i])]
  }
}

sample_normal <- TCGAquery_SampleTypes(colnames(gem_norm[["counts"]]), c("NT"))
exp_matrix_normal <- exp_matrix[,sample_normal]
colnames(exp_matrix_normal) <- substr(colnames(exp_matrix_normal),1,12)
exp_matrix_normal <- t(exp_matrix_normal)
exp_matrix_normal <- exp_matrix_normal[sort(rownames(exp_matrix_normal)),]
write.table(exp_matrix_normal, file= paste0(cancer, "_exp_matrix_normal.tsv"), 
            row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
exp_matrix_normal_log <- log2(exp_matrix_normal+1.0001)
write.table(exp_matrix_normal_log, file= paste0(cancer, "_exp_matrix_normal_log.tsv"), 
            row.names=TRUE, na="NA",col.names=TRUE, sep="\t")

sample_tumor <- TCGAquery_SampleTypes(colnames(gem_norm[["counts"]]), c("TP"))
exp_matrix_tumor <- exp_matrix[,sample_tumor]
colnames(exp_matrix_tumor) <- substr(colnames(exp_matrix_tumor),1,12)
exp_matrix_tumor <- t(exp_matrix_tumor)
exp_matrix_tumor <- exp_matrix_tumor[sort(rownames(exp_matrix_tumor)),]
write.table(exp_matrix_tumor, file= paste0(cancer, "_exp_matrix_tumor.tsv"), 
            row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
exp_matrix_tumor_log <- log2(exp_matrix_tumor+1.0001)
write.table(exp_matrix_tumor_log, file= paste0(cancer, "_exp_matrix_tumor_log.tsv"), 
            row.names=TRUE, na="NA",col.names=TRUE, sep="\t")

save.image()


# Sign Degs with Gene Symbols
sign_degs_2 <- subset(sign_degs_with_ids, select=c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "logFC", "adj.P.Val"))
#save(sign_degs_2, file = paste0(cancer, "_sign_degs_2.rda"))
write.table(sign_degs_2 , file = paste0(cancer, "_sign_degs_2.tsv"), row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

sign_down <- sign_degs_2[(sign_degs_2$adj.P.Val<0.01 & sign_degs_2$logFC<(-1)),]
write.table(sign_down, file = paste0(cancer, "_sign_down.tsv"), row.names=FALSE, na="NA",col.names=TRUE, sep="\t")
sign_up <- sign_degs_2[(sign_degs_2$adj.P.Val<0.01 & sign_degs_2$logFC>1),]
write.table(sign_up, file = paste0(cancer, "_sign_up.tsv"), row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

sign_degs_3 <- subset(sign_degs_2, select="logFC")
sign_degs_3 <- unlist(sign_degs_3, use.names=FALSE)
names(sign_degs_3) <- sign_degs_2$entrezgene_id
sign_degs_3 <- sort(sign_degs_3, decreasing = TRUE)

save.image()

### Enrichment Analysis

library(clusterProfiler)
library(org.Hs.eg.db)

## KEGG over-representation test
# Input ID type can be 'kegg', 'ncbi-geneid', 'ncbi-proteinid' or 'uniprot'.
EA <- enrichKEGG(gene = sign_degs_2$entrezgene_id,
                 keyType = 'ncbi-geneid',
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

EA <- setReadable(EA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA), file = paste0(cancer, "_pathways.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")


EA_up <- enrichKEGG(gene = sign_up$entrezgene_id,
                    keyType = 'ncbi-geneid',
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

EA_up <- setReadable(EA_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA_up), file = paste0(cancer, "_pathways_up.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")


EA_down <- enrichKEGG(gene = sign_down$entrezgene_id,
                      keyType = 'ncbi-geneid',
                      organism = 'hsa',
                      pvalueCutoff = 0.05)

EA_down <- setReadable(EA_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA_down), file = paste0(cancer, "_pathways_down.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

save.image()
}

print(paste0("########## ", cancer, " ", cohort, " DEA ANALYSIS COMPLETED", " ##########"))

save.image()
q(save="yes")

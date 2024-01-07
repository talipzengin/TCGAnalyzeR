### SNV Analysis of TCGA------------------------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(maftools)
library(readr)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least one arguments (TCGA cancer type) must be supplied.n", call.=FALSE)
}

cancer <- args[1]
project <- paste0("TCGA-", cancer)
project_ <- gsub("-","_", project)

data_dir <- "/media/data02/GDC_data/TCGA/"

wd <- paste0(project, "/SNV")
dir.create(wd)
setwd(wd)

pd <- dirname(getwd())
ppd <- dirname(pd)


load(file = paste0(data_dir, project, "/", project_, "_Simple_Nucleotide_Variation_Masked_Somatic_Mutation", ".rda"))


# These codes below determine the sample barcodes of patients in mutation data.
cases_maf <- maf@clinical.data[["Tumor_Sample_Barcode"]]
cases_maf <- as.character(unlist(cases_maf))
cases_maf <- sort(cases_maf)

# Filter TCGA Replicate Samples (cases_maf -> uniq_tsb_maf)
source(paste0(ppd,"/tcga_replicateFilter.R"))
write.table(uniq_tsb_maf, file = paste0(cancer, "_cases_snv.tsv"), 
            row.names=FALSE, na="NA", col.names=FALSE, sep="\t")

# subsetMaf command from maftools package, extracts the mutation data of indicated barcodes, 
maf = subsetMaf(maf, tsb = uniq_tsb_maf, mafObj = TRUE)
#save(maf, file=paste0(cancer, "_maf.RData"))
write.mafSummary(maf = maf, basename = cancer)

save.image()

# The command below summarizes the number, mean and median of mutation types in mutation data as a table.
write.table(maf@summary, file = paste0(cancer, "_maf_summary.tsv"), row.names=FALSE, na="NA", col.names=TRUE, sep="\t")

# The command below summarizes the number of mutation types in each sample as a table.
write.table(maf@variant.type.summary, file = paste0(cancer, "_maf_variant.type.summary.tsv"), row.names=FALSE, na="NA", col.names=TRUE, sep="\t")

## COSMIC Reference Files ---

# First download mutation reference file (CosmicMutantExport.tsv.gz) under "COSMIC Mutation Data" 
#section and list of all cancer census genes (cancer_gene_census.csv) under "Cancer Gene Census" 
#section file from the COSMIC database as a registered user.

# Open these files (cosmic_file specifies the path and the name of the downloaded file)
#cosmic_file <- paste0(ppd, "/CosmicMutantExport.tsv")
#cosmic <- read.table(cosmic_file, header=TRUE, sep='\t', row.names=NULL, stringsAsFactors=FALSE, quote='')
#cosmic <- cosmic %>% select(Gene.name, Gene.CDS.length, Sample.name, Mutation.Description, Mutation.AA, 
#                            Mutation.CDS, Mutation.genome.position) %>% rename(Mutation.GRCh37.genome.position = Mutation.genome.position)


## Preparation of Reference Files for SomInaClust 

#library(SomInaClust)

# This command prepares the reference files which will be used during the determination of significant driver mutations.
#SomInaClust_ref(maf=cosmic,
#                database='COSMIC',
#                convert_cosmic_to_tcga=TRUE,
#                convert_genenames_to_HGNC=TRUE,
#                calculate_CDS=TRUE,
#                convert_cosmic_tcga_input=converttable_mutation_cosmic_tcga,
#                convert_genenames_HGNC_input=converttable_genenames_HGNC,
#                n_cores = 4,
#                filename="SomInaClust_ref.txt")

# This commands to load ref data 
#ref_CDSlength <- read.table("SomInaClust_ref_CDSlength.txt", sep="\t",header=TRUE,
#                            row.names=1,stringsAsFactors = FALSE)
#ref_clusters <- read.table("SomInaClust_ref_clusters.txt", sep="\t",header=TRUE,
#                           stringsAsFactors = FALSE)
#ref_corr_factors <- read.table("SomInaClust_ref_corr_factors.txt", sep="\t",header=TRUE,
#                               row.names=1,stringsAsFactors = FALSE)

# This command to load cancer gene census file
#cancer_gene_census <- read.csv(paste0(ppd, "/cancer_gene_census.csv"), row.names=1)
#cancer_gene_census <- cancer_gene_census %>% select(Entrez.GeneId, Molecular.Genetics)

# Save the ref files to use later and load for another analysis in order to avoid loosing time by preparing ref files again and again
#save(ref_CDSlength, ref_clusters, ref_corr_factors, cancer_gene_census, file="SomInaClust_ref.RData")


## Driver Mutation Analysis by SomInaClust

library(SomInaClust)

# Load cosmic GRCh38 data (to be updated if necessary)
load(paste0(ppd, "/SomInaClust_ref/SomInaClust_ref.rda"))

# Add colname with CDS position
CDS_pos <- maf@data[,"CDS_position"] %>% unlist
CDS_pos <- gsub("-.*|/.*","", CDS_pos)
maf@data[,"CDS_Position"] <- CDS_pos

# Add colname with c_position_WU to use calculate_CDS = TRUE option of SomInaClust_det command
c_position_WU <- paste0('c.', CDS_pos)
maf@data[,"c_position_WU"] <- c_position_WU

# This command
SomInaClust_result <- SomInaClust_det(maf = maf@data, 
                                      database = "TCGA",
                                      save_output = TRUE,
                                      return_results = TRUE,
                                      convert_cosmic_to_tcga = FALSE,
                                      convert_genenames_to_HGNC = FALSE,
                                      calculate_CDS = TRUE,
                                      save_new_maf = TRUE,
                                      CDS_length = ref_CDSlength,
                                      cluster_matrix_input = ref_clusters,
                                      corr_factor_input = ref_corr_factors,
                                      filename = paste0(cancer, "_SomInaClust_results.txt"),
                                      n_cores = 4,
                                      create_pyramidplot = TRUE,
                                      pyramidplot_n_genes = 40,
                                      create_summary_table = FALSE,
                                      CGC = cancer_gene_census,
                                      colname_CGC_class = "Molecular.Genetics")

save.image()

####################

SomInaClust_result <- as.data.frame(SomInaClust_result)


SomInaClust_summary <- function(
    SomInaClust_result,
    save_output=TRUE,
    filename="SomInaClust_results_summary_table.txt",
    save_location=getwd(),
    CGC=NA,
    colname_CGC_class="Molecular.Genetics")
  {
  filename<-paste(save_location,"/",filename,sep="")
  SomInaClust_sel<-SomInaClust_result[SomInaClust_result[,"DG_q"]<=0.05,c("n_mut","DG_q","OG_score","TSG_score")]
  SomInaClust_sel[,"DG_q"]<-signif(SomInaClust_sel[,"DG_q"],3)
  SomInaClust_sel[,"OG_score"]<-round(100*SomInaClust_sel[,"OG_score"],1)
  SomInaClust_sel[,"TSG_score"]<-round(100*SomInaClust_sel[,"TSG_score"],1)
  Class<-NULL
  CGC_class<-NULL
  for(j in 1:nrow(SomInaClust_sel)){
    if(SomInaClust_sel[j,"TSG_score"]>=20) Class<-append(Class,"TSG")
    else if(SomInaClust_sel[j,"OG_score"]>=20) Class<-append(Class,"OG")
    else Class<-append(Class,NA)
    if(!is.na(CGC)&&rownames(SomInaClust_sel)[j]%in%rownames(CGC)){
      if(CGC[rownames(SomInaClust_sel)[j],colname_CGC_class]=="Dom") CGC_class<-append(CGC_class,"Dom")
      else if (CGC[rownames(SomInaClust_sel)[j],colname_CGC_class]=="Rec") CGC_class<-append(CGC_class,"Rec")
      else CGC_class<-append(CGC_class,NA)
      }
    else CGC_class<-append(CGC_class,NA)
  }
  SomInaClust_sel<-cbind(rownames(SomInaClust_sel),SomInaClust_sel,Class,CGC_class)
  colnames(SomInaClust_sel)<-c("Gene","# Mutations","qDG","OG Score","TSG Score","Classification","CGC")
  if(save_output==TRUE) write.table(SomInaClust_sel,filename,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  SomInaClust_sel
}


SomInaClust_results_summary <- SomInaClust_summary(
  SomInaClust_result,
  save_output=TRUE,
  filename=paste0(cancer, "_SomInaClust_results_summary_table.txt"),
  save_location=getwd(),
  CGC = cancer_gene_census,
  colname_CGC_class = "Molecular.Genetics")

###########################################

### Enrichment Analysis

library(clusterProfiler)
library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
library(readr)

## KEGG over-representation test
EA_total <- enrichKEGG(gene = unique(maf@data$Entrez_Gene_Id), # Input ID type can be 'kegg', 'ncbi-geneid', 'ncbi-proteinid' or 'uniprot'.
                 keyType = 'ncbi-geneid',
                 organism = 'hsa',
                 pvalueCutoff = 0.05)


EA_total <- setReadable(EA_total, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA_total), file = paste0(cancer, "_gene_pathways.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

save.image()

###############

if (file.exists(paste0(cancer, "_SomInaClust_results_summary_table.txt"))){
results_summary <- read_delim(paste0(cancer, "_SomInaClust_results_summary_table.txt"), "\t", 
                              escape_double = FALSE, trim_ws = TRUE)
results_summary <- results_summary[!results_summary$Gene=="TTN",]
driver_genes <- results_summary$Gene

driver_genes_ids <- unique(maf@data$Entrez_Gene_Id[which(maf@data$Hugo_Symbol %in% driver_genes)])

## KEGG over-representation test
EA <- enrichKEGG(gene = driver_genes_ids,
                 keyType = 'ncbi-geneid',
                 organism = 'hsa',
                 pvalueCutoff = 0.05)


EA <- setReadable(EA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(as.data.frame(EA), file = paste0(cancer, "_driver_gene_pathways.tsv"),row.names=FALSE, na="NA",col.names=TRUE, sep="\t")

save.image()
}

## Subset OncoPlot for Driver Genes
if (exists("driver_genes")){
maf_driver = subsetMaf(maf, genes = driver_genes, mafObj = TRUE)
#save(maf_driver, file=paste0(cancer, "_maf_driver_genes.RData"))
write.mafSummary(maf = maf_driver, basename = paste0(cancer, "_driver_genes"))

save.image()
}


print(paste0("########## ", cancer, " SNV ANALYSIS COMPLETED", " ##########"))

save.image()
q(save="yes")

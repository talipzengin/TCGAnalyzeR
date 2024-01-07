library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(maftools)
library(readr)
library(tidyverse)

wd <- paste0("SomInaClust_ref")
setwd(wd)

cosmic_dir <- "COSMIC_GRCh38_v95/"

## COSMIC Reference Files ---

#First download mutation reference file (CosmicMutantExport.tsv.gz) under "COSMIC Mutation Data" 
#section and list of all cancer census genes (cancer_gene_census.csv) under "Cancer Gene Census" 
#section file from the COSMIC database as a registered user.

#Open these files (cosmic_file specifies the path and the name of the downloaded file)
cosmic_file <- paste0(cosmic_dir, "CosmicMutantExport.tsv")
cosmic <- read.table(cosmic_file, header=TRUE, sep='\t', row.names=NULL, stringsAsFactors=FALSE, quote='')
cosmic <- cosmic %>% select(Gene.name, Gene.CDS.length, Sample.name, Mutation.Description, Mutation.AA, 
                            Mutation.CDS, Mutation.genome.position) %>% rename(Mutation.GRCh37.genome.position = Mutation.genome.position)

## Driver Mutation Analysis by SomInaClust 

library(SomInaClust)

# This command prepares the reference files which will be used during the determination of significant driver mutations.
SomInaClust_ref(maf=cosmic,
                database='COSMIC',
                convert_cosmic_to_tcga=TRUE,
                convert_genenames_to_HGNC=TRUE,
                calculate_CDS=TRUE,
                convert_cosmic_tcga_input=converttable_mutation_cosmic_tcga,
                convert_genenames_HGNC_input=converttable_genenames_HGNC,
                n_cores = 4,
                filename="SomInaClust_ref.txt")

# This commands to load ref data 
ref_CDSlength <- read.table("SomInaClust_ref_CDSlength.txt", sep="\t",header=TRUE,
                            row.names=1,stringsAsFactors = FALSE)
ref_clusters <- read.table("SomInaClust_ref_clusters.txt", sep="\t",header=TRUE,
                           stringsAsFactors = FALSE)
ref_corr_factors <- read.table("SomInaClust_ref_corr_factors.txt", sep="\t",header=TRUE,
                               row.names=1,stringsAsFactors = FALSE)

# This command to load cancer gene census file
cancer_gene_census <- read.csv(paste0(cosmic_dir, "cancer_gene_census.csv"), row.names=1)
cancer_gene_census <- cancer_gene_census %>% select(Entrez.GeneId, Molecular.Genetics)

# Save the ref files to use later and load for another analysis in order to avoid loosing time by preparing ref files again and again
save(ref_CDSlength, ref_clusters, ref_corr_factors, cancer_gene_census, file="SomInaClust_ref.rda")

print(paste0("########### ", "SNV_ref HAS BEEN PREPARED", " ###########"))

save.image()
q(save="yes")

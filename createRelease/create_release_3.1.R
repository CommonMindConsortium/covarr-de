# Gabriel Hoffman
# January 13, 2021

# Create a file with only DLPFC counts

# source("createRelease/create_release_3.1.R")

library(synapser)
library(dplyr)
library(githubr)

synLogin()

library(data.table)

# Read raw counts from MPP and HBCC
counts.MPP = fread(synGet('syn21867938')$path, data.table=FALSE)
# Get metadata
downloadFile <- function(id){
  fread(synGet(id)$path, data.table = F)
}
downloadFile_version <- function(id , version){
  fread(synGet(id, version = version)$path, data.table = F)
}

ALL_USED_IDs = c()

# Get ancestry vector calculated using gemtools
ANCESTRY_ID = 'syn9922992'
ALL_USED_IDs = c(ALL_USED_IDs, ANCESTRY_ID)
ANCESTRY.HBCC = downloadFile(ANCESTRY_ID) %>% 
  plyr::rename(c(ID = 'SNP_report:Genotyping_Sample_ID'))

ANCESTRY_ID = 'syn2511399'
ALL_USED_IDs[length(ALL_USED_IDs)+1] = ANCESTRY_ID
ANCESTRY.MPP = downloadFile(ANCESTRY_ID) %>% 
  plyr::rename(c('DNA_report..Genotyping.Sample_ID' = 'SNP_report:Genotyping_Sample_ID'))

ANCESTRY = rbind(ANCESTRY.HBCC[,colnames(ANCESTRY.MPP)[-2]], ANCESTRY.MPP[,colnames(ANCESTRY.MPP)[-2]])

# Get genotype ids from synapse. 
GENOTYPE_ID = 'syn16816490'
ALL_USED_IDs = c(ALL_USED_IDs, GENOTYPE_ID)
GENOTYPE = downloadFile(GENOTYPE_ID) %>% 
  dplyr::select(Individual_ID, `SNP_report:Genotyping_Sample_ID`, `SNP_report:Exclude?`) %>% 
  dplyr::inner_join(ANCESTRY) %>% 
  dplyr::filter(is.na(`SNP_report:Exclude?`))

# Get RNASeq QCmetadata
METADATA_QC_DLPFC_ID = 'syn16816488' 
ALL_USED_IDs = c(ALL_USED_IDs, METADATA_QC_DLPFC_ID)
metadata = downloadFile_version(METADATA_QC_DLPFC_ID, version = 19) %>%
  dplyr::select(`Individual_ID`, Institution, Cohort, `Reported_Gender`, Sex, Ethnicity, ageOfDeath, `PMI_(in_hours)`, Dx, 'pH',
                Sample_RNA_ID, one_of('rnaSeq_isolation:RIN',
                                      'rnaSeq_dissection:Brain_Region')) 

# Keep only DLPFC samples
METADATA = metadata %>%
  dplyr::left_join(GENOTYPE) %>% dplyr::rename(Region = `rnaSeq_dissection:Brain_Region`, 
                SampleID = Sample_RNA_ID) %>% dplyr::select(SampleID, Region) %>% dplyr::filter(SampleID %in% colnames(counts.MPP)) %>% 
  		 		 dplyr::filter(Region %in% c('DLPFC')) 

# get counts from only DLPFC
idx = which(colnames(counts.MPP) %in% METADATA$SampleID)
countsFilter = counts.MPP[,c(1:6, idx)]

# write result top file
file = "MSSM.Penn.Pitt_DLPFC.featureCount.tsv"
write.table(countsFilter, file = file,  sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
system(paste("gzip -f", file))


# 
thisFileName <- 'create_release_3.1.R'

# Github link
thisRepo <- getRepo(repository = "CommonMindConsortium/covarr-de", ref="branch", refName='master')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('createRelease/',thisFileName))


nEXP_OBJ = File(paste0(file, ".gz"), 
                name = "Count files for DLFPC from MSSM, Penn, Pitt", 
                parentId = 'syn24172721')
synStore(nEXP_OBJ, used = ALL_USED_IDs, executed = thisFile)



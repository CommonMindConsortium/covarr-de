# Gabriel Hoffman
#
# August 31, 2020
#
# Link FANS4 identifier to Individual_ID for genotype and Sample_RNA_ID for bulk RNA-seq

library(data.table)
library(openxlsx)
library(dplyr)
library(synapser)

synLogin()

# Get meta-data for FANS4 data
###############################

file = synGet('syn22321012')$path
df_featureCounts_gene = fread( file )
geneCounts = df_featureCounts_gene[,7:ncol(df_featureCounts_gene)]

# read metadata
file = synGet('syn22321013')$path
metadata_FANS4 = read.xlsx( file )

# extract information from sample names
info = data.frame( Name = colnames(geneCounts), stringsAsFactors=FALSE)
rownames(info) = info$Name
info$DissectionID = sapply(strsplit( info$Name, '\\.'), function(x) x[1])
info$CellType = sapply(strsplit( info$Name, '\\.'), function(x) x[2])
info$Assay = sapply(strsplit( info$Name, '\\.'), function(x) x[3])

# set color codes for each cell type
colorCodes = c(
  GABA = '#66A61E',
  GLU = '#E6AB02',
  Olig = '#E7298A',
  MgAs = '#7570B3')

# set factor levels to be the same as the order of the colors
# this makes plotting easier
info$CellType = factor(info$CellType,names(colorCodes))

# add information about each sample
info = merge(info, metadata_FANS4, by.x="DissectionID", by.y='Dissection.ID')
rownames(info) = info$Name

# set order to be the same as geneCounts
info = info[match(colnames(geneCounts), info$Name),]

# Intersect with BULK RNA-seq
#############################

# Get bulk RNA-seq counts
file = synGet('syn22045729')$path
geneCountsBulk = data.frame(fread( file ))
rownames(geneCountsBulk) = geneCountsBulk[,1]
geneCountsBulk = data.frame(geneCountsBulk[,-1])

# get link between Individual_ID and Sample_RNA_ID 
METADATA = fread(synGet('syn16816488')$path , data.table=FALSE) %>% 
	dplyr::select(Individual_ID, Sample_RNA_ID)

info2 = merge(info, METADATA, by.x="Individual.ID", by.y="Individual_ID")

# Result linking
# Name (FANS4 sample identifier)
# Individual.ID (identifies a single individual in the genotype file)
# Sample_RNA_ID (identifer for bulk RNA-seq) 
with(info2, unique(cbind(Name, Individual.ID, Sample_RNA_ID)))



#-----------------------------------------------------------------------------------
# preprocessing
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)
rm(list=ls())

library(minfi)
library(GEOquery)
library(limma)
library(openxlsx)

source("./R/MNPprocessIDAT_functions.R")

dir.create("./results")

# get sample annotation from GEO
gse <- getGEO("GSE90496", GSEMatrix=TRUE, getGPL=FALSE)
anno <- pData(gse$GSE90496_series_matrix.txt.gz)

# train just on the ATRTs 
anno <- anno[grep("ATRT",anno$`methylation class:ch1`),]

# read raw data downloaded from GEO and extracted in GSE90496_RAW
filepath <- paste0("GSE90496_RAW/",gsub("_Grn.*","",gsub(".*suppl/","",anno$supplementary_file)))
RGset <- read.metharray(filepath,verbose=TRUE)

# Illumina normalization
message("running normalization ...",Sys.time())
Mset <- MNPpreprocessIllumina(RGset)

# probe filtering
message("probe filtering ...",Sys.time())
amb.filter <- read.table("./filter/amb_3965probes.vh20151030.txt",header=F)
epic.filter <- read.table("./filter/epicV1B2_32260probes.vh20160325.txt",header=F)
snp.filter <- read.table("./filter/snp_7998probes.vh20151030.txt",header=F)
xy.filter <- read.table("./filter/xy_11551probes.vh20151030.txt",header=F)
rs.filter <- grep("rs",rownames(Mset))
ch.filter <- grep("ch",rownames(Mset))

# filter CpG probes
remove <- unique(c(match(amb.filter[,1], rownames(Mset)),
                   match(epic.filter[,1], rownames(Mset)),
                   match(snp.filter[,1], rownames(Mset)),
                   match(xy.filter[,1], rownames(Mset)),
                   rs.filter,
                   ch.filter))

Mset_filtered <- Mset[-remove,]

save(Mset,anno,file="./results/Mset_filtered.RData")  

rm(Mset)
gc()

#batch adjustment
message("performing batchadjustment ...",Sys.time())

methy <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)
rm(Mset_filtered)
gc()

# get FFPE/Frozen type
ffpe <- anno$`material:ch1`
batch <- ifelse(ffpe == "FFPE", 2, 1)

# remove batch effects by linear model
methy.ba <- 2^removeBatchEffect(log2(methy +1), batch)
unmethy.ba <- 2^removeBatchEffect(log2(unmethy +1), batch)

# extract effects to adjust diagnostic samples
s.frozen <- min(which(batch == 1))
s.ffpe <- min(which(batch == 2))
methy.coef <- unmethy.coef <- list()
methy.coef[["Frozen"]] <- log2(methy.ba[, s.frozen]) - log2(methy[, s.frozen] +1)
methy.coef[["FFPE"]] <- log2(methy.ba[, s.ffpe]) - log2(methy[, s.ffpe] +1)
unmethy.coef[["Frozen"]] <- log2(unmethy.ba[, s.frozen]) - log2(unmethy[, s.frozen] +1)
unmethy.coef[["FFPE"]] <- log2(unmethy.ba[, s.ffpe]) - log2(unmethy[, s.ffpe] +1)

# save batch effects 
save(methy.coef,unmethy.coef,file="./results/ba.coef.RData") 

# recalculate betas, illumina like
betas <- methy.ba / (methy.ba +unmethy.ba +100)
betas <- as.data.frame(t(betas))
save(betas,anno,file="./results/betas.ba.RData")  
message("preprocessing finished ...",Sys.time())


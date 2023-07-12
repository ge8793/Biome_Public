#########################################################################################################
######## PROJECT: Running MR analysis between metabolites and Breast Cancer (from BCAC) with Shin et al. metabolite data 
######## Script: Grace Edmunds adapted from Kaitlin Wade, Performing MR analysis with metabolites as exposure and breast cancer as outcome
######## Date: 12/05/23

rm(list=ls())

### INSTALL PACKAGES
# install.packages("devtools")
# # devtools::install_github("MRCIEU/TwoSampleMR") #to update the package
# # devtools::install_github("MRCIEU/MRInstruments")
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("xlsx")
# install.packages("png")
# install.packages("ggrepel")
# install.packages("ggthemes")
# install.packages("calibrate")
# install.packages("data.table")
# install.packages("stringr")
# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR", force = TRUE)
# devtools::install_github("mrcieu/ieugwasr", force = TRUE)
# install.packages("googleAnalyticsR")
library(stringr)
library(remotes)
library(data.table)
library(calibrate)
library(ggrepel)
library(ggthemes)
library(devtools)
library(TwoSampleMR)
library(ieugwasr)
library(MRInstruments)
library(xlsx)
library(plyr) 
library(dplyr)
library(ggplot2)
library(png)
library(tidyr)
library(forestplot)
library(googleAnalyticsR)
ga_auth(email="myemailhere")

# clearws
rm(list=ls())

setwd("mywdhere")


#########################################
## PREPARING EXPOSURE DATA FROM METABOLITE GWAS FILES ##
#########################################

# We downloaded Supplementary table 5 from the Shin et al. GWAS paper then saved as "shin_summary_stats.csv" and read this in and format

Shin <- as.data.frame(fread("./data/shin_summary_stats.csv", skip = 2, select = c(1:10)))
head(Shin)
names(Shin) <- c("Locus", "SNP", "Biochemical", "Ratio", "Effect/Other", "NA?", "n", "eaf", "beta/se", "Pval")
head(Shin)
Shin <- Shin[,c(1:5,7:10)]
Shin <- Shin[c(3:326),]

#Split the columns that are shared
newcols <- str_split_fixed(Shin$`Effect/Other`, "/", 2)
newcols1 <- as.data.frame(newcols)
Shin$effect_allele <- newcols1[,1]
Shin$other_allele <- newcols[,2]

newcols <- str_split_fixed(Shin$`beta/se`, " ", 2)
newcols1 <- as.data.frame(newcols)
Shin$b <- newcols1[,1]
Shin$se <- newcols[,2]
#get rid of brackets by getting the number of characters in the se variable then limiting those you want
nchar(Shin$se)
Shin$se <- as.numeric(substr(Shin$se, start = nchar(Shin$se)-5, stop = nchar(Shin$se)-1))
Shin$se <- str_replace(Shin$se, " ", "")
Shin$n <- str_replace(Shin$n, ",", "")
Shin$n <- as.numeric(Shin$n)
Shin$b <- as.numeric(Shin$b)
Shin$se <- as.numeric(Shin$se)
Shin$Pval <- as.numeric(Shin$Pval)
Shin$eaf <- as.numeric(Shin$eaf)
Shin$SNP <- as.character(Shin$SNP)
Shin$effect_allele <- as.character(Shin$effect_allele)
Shin$other_allele <- as.character(Shin$other_allele)

write.table(Shin, "./data/Shin_Supp5_Formatted.txt", sep = '\t', row.names = F, col.names = T, quote = F)

#Or, if you have done the above before, just read in the exposure data
exposure_dat <- read_exposure_data("./data/Shin_Supp5_Formatted.txt", clump = FALSE, sep="\t", snp_col = "SNP", beta_col = "b", se_col = "se", pval_col = "Pval", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", samplesize_col = "n", phenotype_col = "Biochemical")


#################### PREPARE OUTCOME DATA #####################
#Read in BCAC 20202 Subset BC data and prep it for use =======

#We downloaded the full BCAC outcome data from https://bcac.ccge.medschl.cam.ac.uk/bcacdata/icogs-complete-summary-results/
#In this analysis we are looking at subsetted BC so use the subset dataset saved as 'icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt'

BCAC <- as.data.frame(fread("./data/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt", 
                            select = c(1,2,7,8,9,21,26,27,28,40,41,42:53)))
BCAC2 <- BCAC    
BCAC <- BCAC2
BCACtrim <- BCAC

colnames(BCACtrim) <- c("variant", "snp.icogs", "eaf.controls.icogs", "eaf.cases.iCOGs", "R2.iCOGs", "snp.onco", "eaf.controls.onco", 
                        "eaf.cases.onco", "R2.onco", "effect.allele.outcome", "ref.allele.outcome", "beta.luminalA", "sd.luminalA", "beta.luminalB", "sd.luminalB", "beta.luminalBHER2neg", "sd.luminalBHER2neg",  
                        "beta.HER2enriched", "sd.HER2enriched", "beta.TNBC", "sd.TNBC", "FTOP.pval", "MTOP.pval")
head(BCACtrim)


#Add sample size as described in 2020 paper methods in supplementary table 4 on the line reading subtype analysis
#We contacted the authors to verify this sample size owing to the post hoc analysis used in the meta-GWAS
BCACtrim$ncase <- as.numeric("91562")
BCACtrim$ncontrol <- as.numeric("94407")
BCACtrim$samplesize.outcome <- as.numeric("185969")

#split the snp column to get just the rsid out of it 
BCACtrim[,27:30] <- str_split_fixed(BCACtrim$snp.icogs, ":", 4)
#preview the first ten rows
BCACtrim [1:10,]
#rename the rsid column
colnames(BCACtrim)[27] <- "rsid"
#preview the first ten rows
BCACtrim [1:10,]
#delete the last 3 columns
BCACtrim <- BCACtrim[,1:27]
#preview the first ten rows
BCACtrim [1:10,]
colnames(BCACtrim)
BCAC <- BCACtrim
#trimfurther
BCACtrim <- BCACtrim[,c(1:19, 24:27)]
colnames(BCACtrim)

#Add eaf column using Ryan formula 
class(BCACtrim)
BCACtrim$eaf <- NA
head(BCACtrim)
BCACtrim$eaf <- (((BCACtrim$eaf.cases.iCOGs * 34783) + (BCACtrim$eaf.controls.icogs * 37628) + (BCACtrim$eaf.cases.onco * 71708) + (BCACtrim$eaf.controls.onco * 56779)) / 185969)
BCACtrim [1:10,]

#Make a copy of the data and trim it further
BCAC <- BCACtrim
colnames(BCACtrim)
BCACtrim <- BCAC[,c(5,10:24)]
head(BCACtrim)
colnames(BCACtrim)

#Rename
colnames(BCACtrim) <- c("r2", "effect.allele.outcome", "ref.allele.outcome", "beta.luminalA", "sd.luminalA", "beta.luminalB", "sd.luminalB", "beta.HER2neg", 
"sd.HER2neg","beta.HER2pos", "sd.HER2pos", "ncase.outcome", "ncontrol.outcomee","samplesize.outcome", "SNP", "eaf.outcome")   

#output this new dataframe as a csv file
write.table(BCACtrim, "./data/BCAC_subtypes_2020.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t') 

#Read back in trimmed data when needed 
BCACtrim <- as.data.frame(fread("./data/BCAC_subtypes_2020.txt"))

#Read in metabolite exposure data
exposure_dat <- read_exposure_data("./data/Shin_Supp5_Formatted.txt", clump = FALSE, sep="\t", snp_col = "SNP", beta_col = "b", se_col = "se", pval_col = "Pval", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", samplesize_col = "n", phenotype_col = "Biochemical")

metab_snps <- as.character(exposure_dat$SNP)

cancer <- filter(BCACtrim, BCACtrim$SNP %in% metab_snps)
table(cancer$SNP)

write.table(cancer, "./data/Outcome_subtype_bcac_metabsnps.txt", row.names = F, col.names = T, quote = F, sep = '\t')

#If you would like to use exactly this outcome data, it is available on github at https://github.com/ge8793/Biome_Public/blob/main/outcome_subtype_bcac_metabsnps.txt

#Read in outcome data

cancer <- as.data.frame(fread("./data/Outcome_subtype_bcac_metabsnps.txt"))

#make a copy of the data 
cancer2 <- cancer
#preview it
head(cancer)
colnames(cancer)
#select only the columns we want
cancer_luminalA <- cancer [,c(1:3,4,5,12:16)]
cancer_luminalA$outcome <- "Luminal A Breast Cancer"
colnames(cancer_luminalA)[4] <- "beta.outcome"
colnames(cancer_luminalA)[5] <- "se.outcome"
cancer_luminalB <- cancer [,c(1:3,6,7,12:16)]
cancer_luminalB$outcome <- "Luminal B Breast Cancer"
colnames(cancer_luminalB)[4] <- "beta.outcome"
colnames(cancer_luminalB)[5] <- "se.outcome"
cancer_luminalBHER2neg <- cancer [,c(1:3,8,9,12:16)]
cancer_luminalBHER2neg$outcome <- "HER2 negative Breast Cancer"
colnames(cancer_luminalBHER2neg)[4] <- "beta.outcome"
colnames(cancer_luminalBHER2neg)[5] <- "se.outcome"
cancer_HER2enriched <- cancer [,c(1:3,10,11,12:16)]
cancer_HER2enriched$outcome <- "HER2 positive Breast Cancer"
colnames(cancer_HER2enriched)[4] <- "beta.outcome"
colnames(cancer_HER2enriched)[5] <- "se.outcome"


cancer_list <- list(cancer_luminalA, cancer_luminalB, cancer_luminalBHER2neg, cancer_HER2enriched)

# for(i in 1:length(cancer_list)){
# name <- cancer_list[[i]]$outcome[1]
# filepath = paste("./analysis/Metabolites/230703", name, ".txt")
# write.table(cancer_list[[i]], filepath, col.names = T, row.names = F, quote = T, sep = "\t")
# }

# Set the loop up for main analyses
# Making lists to input results within
results_list <- list()
het_results_list <- list()
pleio_results_list <- list()
dat_list <- list ()

colnames(cancer)

# Run the MR analyses for Luminal A

  for(j in 1:length(cancer_list)){
    name <- cancer_list[[j]]$outcome[1]
    filepath = paste("./analysis/Metabolites/160703", name, ".txt")
    outcome_dat <- read_outcome_data(filepath, sep="\t", snps = exposure_dat$SNP,
                                     snp_col = "SNP", beta_col = "beta.outcome", se_col = "se.outcome", 
                                     effect_allele_col = "effect.allele.outcome", other_allele_col = "ref.allele.outcome", 
                                     ncase_col = "ncase.outcome", ncontrol_col = "ncontrol.outcome", eaf_col = "eaf.outcome",
    )
    
    dat <- harmonise_data(exposure_dat, outcome_dat, action =2)
    mr_results <- mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
    mr_results$outcome <- paste0(name)
    # Outputting each MR result into results_list
    a <- cbind.data.frame(mr_results$exposure, mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)
    k <- paste("res_","_",name,sep="")
    l <- paste("dat_", "_", name, sep = "")
    results_list[[paste(name,sep="")]] <- assign(k,a)
    dat_list[[paste(name,sep = "")]] <- assign(l,dat)
    rm(a,k,l,dat,name)
    rm(a,k, name)
    
  }

# Merge all results files together and input into an Excel spreadsheet
all_results <- do.call(rbind, results_list)
rownames(all_results) <- NULL
names(all_results)
names(all_results) <- c("exposure", "outcome", "nsnp", "method", "b", "se", "pvalue")
head(all_results)
mr_results <- all_results
mr_results$OR <- exp(mr_results$b)
mr_results$LCL <- exp((mr_results$b)-1.96*(mr_results$se))
mr_results$UCL <- exp((mr_results$b)+1.96*(mr_results$se))
write.table(mr_results, "./analysis/Metabolites/Subtype_Results/Subtype_Results.txt", row.names=F, col.names = T, quote = F, sep = "\t")

F_list <- list()

for(i in 1:length(dat_list))
{
  dat <- dat_list[[i]] #harmonization
  eaf <- dat$eaf.exposure     
  b <- dat$beta.exposure
  se <- dat$se.exposure
  p <- dat$pval.exposure
  snp<- dat$SNP
  N <- dat$samplesize.exposure
  samplesize <- dat$samplesize.exposure
  k<-1
  dat$VG <- 2*(b^2)*eaf*(1-eaf)
  dat$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
  dat$R2 <- dat$VG/dat$PV
  dat$Fstat = dat$R2*(N-1-k)/((1-dat$R2)*k)
  
  name <- paste0(dat_list[[i]]$exposure[1])
  F_list[[paste(i, "_", dat_list[[i]]$exposure[1], sep = "")]] <- assign(name,dat)
  filename <- paste0("./analysis/", name, "_","harmonised_data",".txt")
  write.table(dat, file = filename, col.names = T, row.names = F, quote = FALSE, sep = "\t")
  rm(name)
}

F_statistics <- rbindlist(F_list)
F_statistics_new <- cbind(F_statistics$exposure, F_statistics$R2, F_statistics$Fstat)
names <- c("Exposure", "R2", "F_statistic")
colnames(F_statistics_new) <- names

write.table(F_statistics, "./analysis/metabolites/F_R2_Power_metabolites_SubtypeBC.txt", col.names = T, row.names = F, quote = FALSE, sep = '\t')

Q_list <- list()

for(i in 1:length(dat_list))if(nrow(dat_list[[i]])>1){
  Q <- mr_heterogeneity(dat_list[[i]])
  Q_baso <- as.data.frame(Q)
  Q$outcome <- paste0("TNBC")
  name <- paste0(dat_list[[i]]$exposure[1])
  Q_list[[paste(i, "_", dat_list[[i]]$exposure[1], sep = "")]] <- assign(name,Q_baso)
  rm(Q_baso, name)
}

Q_statistics <- rbindlist(Q_list)

write.table(Q_statistics, "./analysis/metabolites/Q_metabolites_SubtypeBC.txt", col.names = T, row.names = F, quote = FALSE, sep = '\t')

#Get the MR Egger intercept for carnitine in each case
dat <- dat_list[[1]]
het <- mr_pleiotropy_test(dat)

dat <- dat_list[[2]]
het2 <- mr_pleiotropy_test(dat)

dat <- dat_list[[3]]
het3 <- mr_pleiotropy_test(dat)

dat <- dat_list[[4]]
het4 <- mr_pleiotropy_test(dat)

dat <- dat_list[[5]]
het5 <- mr_pleiotropy_test(dat)

het_all <- rbind(het, het2, het3, het4, het5)

write.table(het_all, "./analysis/metabolites/EggerIntercept_metabolites_SubtypeBC.txt", col.names = T, row.names = F, quote = FALSE, sep = '\t')

#check unique SNPS in exposure data
SNPS <- as.data.frame(unique(exposure_dat$SNP)) #218 SNPs

#check unique SNPs in outcome data 
SNPS <- as.data.frame(unique(cancer$SNP)) #204 SNPs  

#Make a list of the SNPs that are missing and manually look them up in LDlink
proxy_list <- filter(exposure_dat, !exposure_dat$SNP %in% cancer$SNP) #45 snps
sanity <- unique(proxy_list$SNP) # 14 exposures 
#harmonise
dat <- harmonise_data(exposure_data, outcome_data, action =2)
check <- unique(dat$SNP) #204 SNPs 
check <- unique(dat$exposure) #209 exposures 

#limit my results to IVW and Wald only
Res <- all_results %>%
  filter(all_results$method %in% c("Wald ratio", "Inverse variance weighted"))

#See how many exposures for Bonferroni
check <- unique(Res$exposure)


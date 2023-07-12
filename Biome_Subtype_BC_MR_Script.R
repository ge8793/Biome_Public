#########################################################################################################
######## PROJECT: Running bi-directional MR analysis between FGFP abundances and Breast Cancer (from BCAC) with new FGFP data
######## Script: Grace Edmunds adapted from Kaitlin Wade, Performing MR analysis with FGFP abundances as exposure and subtype specific breast cancer as outcome
######## Date: 24/06/20

rm(list=ls())

### INSTALL PACKAGES 
# install.packages("devtools")
# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("MRCIEU/MRInstruments")
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("xlsx")
#install.packages("png")
#install.packages("ggrepel")
#install.packages("ggthemes")
#install.packages("calibrate")
#install.packages("data.table")
#install.packages("stringr")
library(stringr)
library(data.table)
library(calibrate)
library(ggrepel)
library(ggthemes)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(xlsx)
library(plyr) 
library(dplyr)
library(ggplot2)
library(png)
library(tidyr)
library(forestplot)

###Check the IEU tools are working as expected 

##check the package is installed 
# devtools::install_github("mrcieu/ieugwasr")

##check you have an access token
ieugwasr::check_access_token()

#if not working, reset email in googleauth
# gar_auth(email = "myemail")

##get an access token if you do not have one
# ieugwasr::get_access_token

##Use this token to see some available outcome data, to confirm it is working 
ao <- available_outcomes()

### Prepare environment

## clearworkspace
rm(list=ls())

##set working directory 
# setwd("<<your directory here>>")

#########################################
## PREPARING GUT MICROBIOME GWAS FILES ##
#########################################
#NB: The bugs specified are those that are presented in Table 1 of the paper as those with associated SNPs.
#NB: Dave created the files that provide the genome-wide results and those considered 'significant' (n=26 across the 19 phenotypes)

#Assign names to run analyses with
bugs <- c("C_Gammaproteobacteria_RNT","F_Sutterellaceae_HB","G_Bifidobacterium_RNT","G_Butyricicoccus_RNT","G_Dialister_HB","G_Parabacteroides_RNT",
          "G_Ruminococcus_HB","G_unclassified_F_Erysipelotrichaceae_HB","G_unclassified_F_Porphyromonadaceae_RNT","G_unclassified_O_Bacteroidales_HB","G_unclassified_P_Firmicutes_HB",
          "G_unclassified_P_Firmicutes_RNT","G_Veillonella_HB","G_Coprococcus_HB")

# Making lists to input results within
results_list <- list()

# Pull out all SNPs that are significant across all exposures
for(i in bugs){
  exposure <- paste("/Users/ge8793/OneDrive - University of Bristol/001_PROJECTS/Biome/data/Bugs/",i,".txt",sep="")
  exposure_dat <- read_exposure_data(exposure,sep="\t", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "P_value", eaf_col = "eaf", effect_allele_col = "allele_B", other_allele_col = "allele_A", samplesize_col = "n")
  exposure_dat$exposure <- i
  assign(paste("exp_meta_",i,sep=""),exposure_dat)
  results_list[[i]] <- assign(i,exposure_dat)
}
all_exposures <- do.call(rbind, results_list)
all_SNPs <- as.character(all_exposures$SNP)
all_SNPs <- as.data.frame(all_SNPs)
colnames(all_SNPs) <- "SNP"

# export this snplist as a table
# write.table(all_SNPs, "./data/bugs_snp_list.txt", row.names= F, col.names = T, quote = F, sep = '\t') 

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

#Add eaf column 
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

cancer <- filter(BCACtrim, BCACtrim$SNP %in% all_SNPs$SNP)
table(cancer$SNP)

write.table(cancer, "./data/outcome_subtype_BC_biome.txt", row.names = F, col.names = T, quote = F, sep = '\t')

#If you would like to use exactly this outcome data, it is available on github at https://github.com/ge8793/Biome_Public/blob/main/outcome_subtype_bcac_metabsnps.txt
#NB: GWAS of breast cancer executed by zhang et al. 2020 using data from BCAC consortium 

# I pulled the bug snps from the overall bcac data using the script 'Biome_from_BCAC.R' 
cancer <- as.data.frame(fread("./data/outcome_subtype_BC_biome.txt"))
#make a copy of the data in case anything goes wrong
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
# 
for(i in 1:length(cancer_list)){
  name <- cancer_list[[i]]$outcome[1]
  filepath = paste("./analysis/", name, ".txt")
  write.table(cancer_list[[i]], filepath, col.names = T, row.names = F, quote = T, sep = "\t")
}

# #There are 3 missing snps, I looked for a proxy using LDLink website, none of the proxies were in the subtype data
#Therefore we cannot instrument
#F_Sutterellaceae_HB
#G_Ruminococcus_HB
#G_Coprococcus_HB


####################################################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND BC (OUTCOME) USING ALL META-SUPPORTED SNPS ##
####################################################################################################
#Assign names to run analyses with
bugs <- c("C_Gammaproteobacteria_RNT","G_Bifidobacterium_RNT","G_Butyricicoccus_RNT","G_Dialister_HB","G_Parabacteroides_RNT"
          ,"G_unclassified_F_Erysipelotrichaceae_HB","G_unclassified_F_Porphyromonadaceae_RNT","G_unclassified_O_Bacteroidales_HB","G_unclassified_P_Firmicutes_HB",
          "G_unclassified_P_Firmicutes_RNT","G_Veillonella_HB")

# Set the loop up for main analyses
# Making lists to input results within
results_list <- list()
het_results_list <- list()
pleio_results_list <- list()
dat_list<- list()

# Run the MR analyses 
for(i in bugs){
  for(j in 1:length(cancer_list)){
  exposure <- paste("/Users/ge8793/OneDrive - University of Bristol/001_PROJECTS/Biome/data/Bugs/",i,".txt",sep="")
  exposure_dat <- read_exposure_data(exposure,sep="\t", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "P_value", eaf_col = "eaf", effect_allele_col = "allele_B", other_allele_col = "allele_A", samplesize_col = "n")
  exposure_dat$exposure <- i
  name <- cancer_list[[j]]$outcome[1]
  filepath = paste("./analysis/", name, ".txt")
  outcome_dat <- read_outcome_data(filepath, sep="\t", snp_col = "SNP", 
                                     beta_col = "beta.outcome", se_col = "se.outcome", pval_col = "pval.outcome", 
                                     eaf_col = "EAF", effect_allele_col = "effect.allele.outcome", 
                                     other_allele_col = "ref.allele.outcome", phenotype_col = "outcome", 
                                     ncase_col = "ncase",
                                     ncontrol_col = "ncontrol",
                                     samplesize_col = "samplesize.outcome",
                                    )

  dat <- harmonise_data(exposure_dat, outcome_dat, action =2)
  filepath = paste0("./analysis/", name, "dat", ".txt")
  write.table(dat, filepath, row.names = F, col.names = T, quote = F, sep = '\t')
  mr_results <- mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
  
  # Outputting each MR result into results_list
  a <- cbind.data.frame(mr_results$exposure, mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)
  k <- paste("res_",i,"_",name,sep="")
  l <- paste("dat_", i, "_", name, sep = "")
  results_list[[paste(i,"_",name,sep="")]] <- assign(k,a)
  dat_list[[paste(i, "_", name,sep = "")]] <- assign(l,dat)
  rm(a,k,l,dat,name)

  }
}

# Merge all results files together and input into an Excel spreadsheet
all_results <- do.call(rbind, results_list)
rownames(all_results) <- NULL
write.table(all_results, "./analysis/mr_results_BCsubtypes2020_biome.txt", row.names=F, col.names = T, quote = F, sep = "\t")

# all_results <- as.data.frame(fread("./analysis/mr_results_BCsubtypes2020_biome.txt"))

###Get the binary F-stats and R2 for this analysis

##Add the prevalence to the binary phenotypes in the dat list from FGFP supplementary data

#C_Gammaproteobacteria_RNT
for(i in 1:5){
dat_list[i] #preview the data
dat_list[[i]]$Ncase.exp <- 1626
dat_list[[i]]$Ncontrol.exp <- 1626
dat_list[[i]]$samplesize.exposure <- 1626
dat_list[[i]]$prevalence <- 1
}

#G_Bifidobacterium_RNT_overall_BC
for(i in 6:10){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 1975
dat_list[[i]]$Ncontrol.exp <- 1975
dat_list[[i]]$samplesize.exposure <- 1975
dat_list[[i]]$prevalence <- 1
}

#G_Butyricicoccus_RNT_overall_BC
for(i in 11:15){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 2223
dat_list[[i]]$Ncontrol.exp <- 2223
dat_list[[i]]$samplesize.exposure <- 2223
dat_list[[i]]$prevalence <- 1
}

#G_Dialister_HB_overall_BC
for(i in 16:20){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 2223*0.50
dat_list[[i]]$Ncontrol.exp <- 2223*0.50
dat_list[[i]]$samplesize.exposure <- 2223
dat_list[[i]]$prevalence <- 0.50
}
#G_Parabacteroides_RNT_overall_BC
for(i in 21:25){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 2223
dat_list[[i]]$Ncontrol.exp <- 2223
dat_list[[i]]$samplesize.exposure <- 2223
dat_list[[i]]$prevalence <- 1
}

#G_unclassified_F_Erysipelotrichaceae_HB_overall_BC
for(i in 26:30){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 2223*0.47
dat_list[[i]]$Ncontrol.exp <- 2223*0.53
dat_list[[i]]$samplesize.exposure <- 2223
dat_list[[i]]$prevalence <- 0.47
}
#G_unclassified_F_Porphyromonadaceae_RNT_overall_BC
for(i in 31:35){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 1420
dat_list[[i]]$Ncontrol.exp <- 1420
dat_list[[i]]$samplesize.exposure <- 1420
dat_list[[i]]$prevalence <- 1
}
#G_unclassified_O_Bacteroidales_HB_overall_BC
for(i in 36:40){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 2223*0.19
dat_list[[i]]$Ncontrol.exp <- 2223*0.81
dat_list[[i]]$samplesize.exposure <- 2223
dat_list[[i]]$prevalence <- 0.19
}
#G_unclassified_P_Firmicutes_HB_overall_BC
for(i in 41:45){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 2223*0.85
dat_list[[i]]$Ncontrol.exp <- 2223*0.15
dat_list[[i]]$samplesize.exposure <- 2223
dat_list[[i]]$prevalence <- 0.85
}
#G_unclassified_P_Firmicutes_RNT_overall_BC
for(i in 46:50){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 1904
dat_list[[i]]$Ncontrol.exp <- 1904
dat_list[[i]]$samplesize.exposure <- 1904
dat_list[[i]]$prevalence <- 1
}
#G_Veillonella_HB_overall_BC
for(i in 51:55){
dat_list[i]
dat_list[[i]]$Ncase.exp <- 2223*0.27
dat_list[[i]]$Ncontrol.exp <- 2223*0.73
dat_list[[i]]$samplesize.exposure <- 2223
dat_list[[i]]$prevalence <- 0.27
}

## Calculate R2 for binary traits (for this we can use the get_r_from_lor function)
res <- all_results
full_dat <- rbindlist(dat_list)
res$exposure <- res$`mr_results$exposure`
res$outcome <- res$`mr_results$outcome`
res2 <- merge(res, full_dat, by = c("exposure", "outcome"))
res2$prevalence <- 1
names(res2)
res2$R <- TwoSampleMR::get_r_from_lor(res2$beta.exposure, res2$eaf.exposure, res2$Ncase.exp, res2$Ncontrol.exp, res2$prevalence, model = "logit", correction = F)
res2$R2 <- res2$R^2
write.table(res2, "./analysis/Subtypes_R2values_binary.txt", row.names = F, col.names = T, quote = F, sep = '\t')

### Get the continuous F-stats and R2 values for this analysis

F_list <- list()

for(i in 1:length(dat_list)){
  dat <- dat_list[[i]] #harmonization
  eaf <- dat$eaf.exposure
  b <- dat$beta.exposure
  se <- dat$se.exposure
  p <- dat$pval.exposure
  snp<- dat$SNP
  N <- (dat$samplesize.exposure)
  samplesize <- dat$samplesize.exposure
  Ncontrols <- dat$Ncontrol.exp
  Ncases <- dat$Ncase.exp
  #individual R2
  k<-1
  n<-dat$samplesize.exposure
  R2 = (2*(b^2)*eaf*(1-eaf))/(2*(b^2)*eaf*(1-eaf))+((se^2)*(2*N)*eaf*(1-eaf))
  dat$R2 <- paste0(R2)
  
  #individual F-stats
    F<-unlist(lapply(1:length(R2),FUN=function(x) R2[x]*(N-1-k)/((1-R2[x])*k)))
    F_stat <- as.character(F)
    dat$F_stat <- F_stat
    rm(F)
  
  #individual power calculation
  sig<-0.05 #alpha
  N<-Ncontrols+Ncases
  ratio<-Ncontrols/Ncases
  hypothesised_b <- 0.1823216 #ln(OR) 
  power<-pnorm(sqrt(N*R2*(ratio/(1+ratio))*(1/(1+ratio)))*hypothesised_b-qnorm(1-sig/2))
  dat$power <- paste0(power)
  dat1 <- as.data.frame(dat)
  name <- paste0(dat_list[[i]]$exposure[1])
  F_list[[paste(i, "_", dat_list[[i]]$exposure[1], sep = "")]] <- assign(name,dat)
  filename <- paste0("./analysis/", dat$exposure[1], "_","harmonised_data",".txt")
  write.table(dat, file = filename, col.names = T, row.names = F, quote = FALSE, sep = "\t")
  rm(dat, name)
}

F_statistics <- rbindlist(F_list)
F_statistics <- cbind(F_statistics$exposure, F_statistics$R2, F_statistics$F_stat, F_statistics$power)
names <- c("Exposure", "R2", "F_statistic", "Power")
colnames(F_statistics) <- names
F_statistics <- unique(F_statistics)

write.table(F_statistics, "./analysis/Continuous_F_R2_Power_MT_SubtypeBC.txt", col.names = T, row.names = F, quote = FALSE, sep = '\t')

#########################################################################################################
######## PROJECT: Running  MR analysis between FGFP abundances and Breast Cancer (from BCAC) with new FGFP data
######## Script: Grace Edmunds, adapted from Kaitlin Wade, Performing MR analysis with FGFP abundances as exposure and breast cancer as outcome
######## Date: 24/06/20


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
ga_auth(email="myemail")

# clearws
rm(list=ls())

getwd()
# setwd("mywd")

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

#############################################################
# PREPARING OUTCOME DATA BY READING in BCAC 20202 CIMBA BC data and prepping it for use =======
##############################################################

#We downloaded the full BCAC outcome data from https://bcac.ccge.medschl.cam.ac.uk/bcacdata/icogs-complete-summary-results/
#In this analysis we are looking at meta-analysed TNBC so use the CIMBA dataset saved as 'CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics 2.txt'

BCAC <- as.data.frame(fread("./data/CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics 2.txt",
                            select = c(2,5,6,7,9,10)))

#see what columns are called
colnames(BCAC)
head(BCAC)
BCACview <- BCAC[1:10,]

#read in a lookup table to convert bugs rsids to chr_posn and then pull them from bcac, this is available on request from the BCAC authors

lookup <- as.data.frame(fread("./data/CIMBA_lookuptable.txt"))

#see what columns are called
colnames(lookup)
lookup[1:10,]

#make the column names into the first row
colnames(lookup) <- lookup[1,]
colnames(lookup) <- c("variant_ID", "SNP")

#pull the biome snps out of the lookup
look_filter <- filter(lookup, lookup$SNP %in% all_SNPs$SNP)

#pull these snps out of cimba
BCAC[1:10,]
cimba_filter <- filter(BCAC, BCAC$MarkerName %in% look_filter$variant_ID)

colnames(cimba_filter)

#Trim the data
BCACtrim <- cimba_filter

colnames(BCACtrim) <- c("variant_ID", "ref.allele.outcome", "effect.allele.outcome", "eaf.outcome", "beta", "se")
head(BCACtrim)

#add the rsids
merge <- merge(BCACtrim, lookup, by = "variant_ID" )

names(merge)
head(merge)

#output this new dataframe as a csv file
write.table(merge, "./data/outcome_CIMBA_2020_Biome.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t') 

#If you want to use this exact outcome data, you can download it from https://github.com/ge8793/Biome_Public/blob/main/outcome_CIMBA_2020_Biome.txt and read in
merge <- as.data.frame(fread("./data/outcome_CIMBA_2020_Biome.txt"))

#add sample size as described in 2020 paper methods 
merge$ncase <- as.numeric("16592") 
merge$ncontrol <- as.numeric("103901") 
merge$samplesize.outcome <- as.numeric("120493") 

names(merge)

outcome_data <- TwoSampleMR::format_data(merge, type = "outcome", snps = all_SNPs$SNP,
                                         snp_col = "SNP", beta_col = "beta", se_col = "se", 
                                         effect_allele_col = "effect.allele.outcome", other_allele_col = "ref.allele.outcome", 
                                         ncase_col = "ncase", ncontrol_col = "ncontrol", eaf_col = "eaf.outcome", samplesize_col = "samplesize.outcome", pval_col = "pval")


cancer <- outcome_data
#make a copy of the data in case anything goes wrong
cancer2 <- cancer
#preview it
head(cancer)
colnames(cancer)

cancer$outcome <- "Meta-analysed Triple Negative Breast Cancer"

##There are 1 missing snps, I looked for a proxy using LDLink website but it was not available  in the CIMBA data so we cannot instrument G.Coprococcus 

####################################################################################################
## RUNNING ANALYSES OF FGFP ABUNDANCES (EXPOSURE) AND BC (OUTCOME) USING ALL META-SUPPORTED SNPS ##
####################################################################################################
#Assign names to run analyses with
bugs <- c("C_Gammaproteobacteria_RNT","F_Sutterellaceae_HB","G_Bifidobacterium_RNT","G_Butyricicoccus_RNT","G_Dialister_HB","G_Parabacteroides_RNT",
          "G_Ruminococcus_HB","G_unclassified_F_Erysipelotrichaceae_HB","G_unclassified_F_Porphyromonadaceae_RNT","G_unclassified_O_Bacteroidales_HB","G_unclassified_P_Firmicutes_HB",
          "G_unclassified_P_Firmicutes_RNT","G_Veillonella_HB")

# Set the loop up for main analyses
# Making lists to input results within
results_list <- list()
het_results_list <- list()
pleio_results_list <- list()
dat_list <- list()

# Run the MR analyses for Luminal A
# Run the MR analyses
for(i in bugs){
  exposure <- paste("./data/Bugs/",i,".txt",sep="")
  exposure_dat <- read_exposure_data(exposure,sep="\t", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "P_value", eaf_col = "eaf", effect_allele_col = "allele_B", other_allele_col = "allele_A", samplesize_col = "n")
  exposure_dat$exposure <- i
  outcome_dat <- read_outcome_data("./data/outcome_CIMBA_2020_Biome.txt", sep="\t", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf.outcome", effect_allele_col = "effect.allele.outcome", other_allele_col = "ref.allele.outcome", phenotype_col = "outcome")
  outcome_dat$samplesize_col <- 120493
  outcome_dat$ncase_col <- 16592 
  outcome_dat$ncontrol_col <- 103901
  dat <- harmonise_data(exposure_dat, outcome_dat, action =2)
  mr_results <- mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
  
  # Outputting each MR result into results_list
  a <- cbind.data.frame(mr_results$exposure, mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)
  name <- paste("res_",i,"_","TNBC_meta",sep="")
  results_list[[paste(i,"_","TNBC_meta",sep="")]] <- assign(name,a)
  dat_list[[paste(i,"_","overall_BC",sep="")]] <- assign(name,dat)
  rm(a,name)
  
}

# Merge all results files together and input into an Excel spreadsheet
all_results <- do.call(rbind, results_list)
rownames(all_results) <- NULL
#write.table(all_results, "./analysis/mr_results_BC_CIMBA2020_biome.txt", row.names=F, col.names = T, quote = F, sep = "\t")

##Add the prevalence to the binary phenotypes in the dat list from FGFP supplementary data

#C_Gammaproteobacteria_RNT
dat_list[1] #preview the data
dat_list[[1]]$Ncase.exp <- 1626
dat_list[[1]]$Ncontrol.exp <- 1626
dat_list[[1]]$samplesize.exposure <- 1626
dat_list[[1]]$prevalence <- 1

#F_Sutterellaceae_HB_overall_BC
dat_list[2]
dat_list[[2]]$Ncase.exp <- 2223*0.93
dat_list[[2]]$Ncontrol.exp <- 2223*0.07
dat_list[[2]]$samplesize.exposure <- 2223
dat_list[[2]]$prevalence <- 0.93

#G_Bifidobacterium_RNT_overall_BC
dat_list[3]
dat_list[[3]]$Ncase.exp <- 1975
dat_list[[3]]$Ncontrol.exp <- 1975
dat_list[[3]]$samplesize.exposure <- 1975
dat_list[[3]]$prevalence <- 1

#G_Butyricicoccus_RNT_overall_BC
dat_list[4]
dat_list[[4]]$Ncase.exp <- 2223
dat_list[[4]]$Ncontrol.exp <- 2223
dat_list[[4]]$samplesize.exposure <- 2223
dat_list[[4]]$prevalence <- 1

#G_Dialister_HB_overall_BC
dat_list[5]
dat_list[[5]]$Ncase.exp <- 2223*0.50
dat_list[[5]]$Ncontrol.exp <- 2223*0.50
dat_list[[5]]$samplesize.exposure <- 2223
dat_list[[5]]$prevalence <- 0.50

#G_Parabacteroides_RNT_overall_BC
dat_list[6]
dat_list[[6]]$Ncase.exp <- 2223
dat_list[[6]]$Ncontrol.exp <- 2223
dat_list[[6]]$samplesize.exposure <- 2223
dat_list[[6]]$prevalence <- 1

#G_Ruminococcus_HB_overall_BC
dat_list[7]
dat_list[[7]]$Ncase.exp <- 2223*0.94
dat_list[[7]]$Ncontrol.exp <- 2223*0.06
dat_list[[7]]$samplesize.exposure <- 2223
dat_list[[7]]$prevalence <- 0.94

#G_unclassified_F_Erysipelotrichaceae_HB_overall_BC
dat_list[8]
dat_list[[8]]$Ncase.exp <- 2223*0.47
dat_list[[8]]$Ncontrol.exp <- 2223*0.53
dat_list[[8]]$samplesize.exposure <- 2223
dat_list[[8]]$prevalence <- 0.47

#G_unclassified_F_Porphyromonadaceae_RNT_overall_BC
dat_list[9]
dat_list[[9]]$Ncase.exp <- 1420
dat_list[[9]]$Ncontrol.exp <- 1420
dat_list[[9]]$samplesize.exposure <- 1420
dat_list[[9]]$prevalence <- 1

#G_unclassified_O_Bacteroidales_HB_overall_BC
dat_list[10]
dat_list[[10]]$Ncase.exp <- 2223*0.19
dat_list[[10]]$Ncontrol.exp <- 2223*0.81
dat_list[[10]]$samplesize.exposure <- 2223
dat_list[[10]]$prevalence <- 0.19

#G_unclassified_P_Firmicutes_HB_overall_BC
dat_list[11]
dat_list[[11]]$Ncase.exp <- 2223*0.85
dat_list[[11]]$Ncontrol.exp <- 2223*0.15
dat_list[[11]]$samplesize.exposure <- 2223
dat_list[[11]]$prevalence <- 0.85

#G_unclassified_P_Firmicutes_RNT_overall_BC
dat_list[12]
dat_list[[12]]$Ncase.exp <- 1904
dat_list[[12]]$Ncontrol.exp <- 1904
dat_list[[12]]$samplesize.exposure <- 1904
dat_list[[12]]$prevalence <- 1

#G_Veillonella_HB_overall_BC
dat_list[13]
dat_list[[13]]$Ncase.exp <- 2223*0.27
dat_list[[13]]$Ncontrol.exp <- 2223*0.73
dat_list[[13]]$samplesize.exposure <- 2223
dat_list[[13]]$prevalence <- 0.27

# Merge all results files together and input into an Excel spreadsheet
all_results <- do.call(rbind, results_list)
rownames(all_results) <- NULL

# write.table(all_results, "./analysis/mr_results_BCoverall2020_biome.txt", row.names=F, col.names = T, quote = F, sep = "\t")

### Get the F-stats and R2 values for this analysis

F_list <- list()

## Calculate R2 for binary traits (for this we can use the get_r_from_lor function)
res <- all_results
full_dat <- rbindlist(dat_list)
res$exposure <- res$`mr_results$exposure`
res2 <- merge(res, full_dat, by = "exposure")
res2$prevalence <- 1
names(res2)
res2$R <- TwoSampleMR::get_r_from_lor(res2$beta.exposure, res2$eaf.exposure, res2$Ncase.exp, res2$Ncontrol.exp, res2$prevalence, model = "logit", correction = F)
res2$R2 <- res2$R^2
R2*(N-1-k)/((1-r2)*k)
res2$F_stat<- (res2$R2*(res2$samplesize.exposure-1))/((1-res2$R2)*1)
# write.table(res2, "./analysis/R2values_binary_CIMBA.txt", row.names = F, col.names = T, quote = F, sep = '\t')
# write.table(res2, "./manuscript/F stats and R2/Binary_CIMBA_R2.txt", row.names = F, col.names = T, quote = F, sep = '\t')

##Calculate R2 for continuous traits 
### For this, we use the equation R2 = 2*(b^2)*eaf*(1-eaf)/(2*(b^2)*eaf*(1-eaf)+(se^2)*(2*N)*eaf*(1-eaf)
### where b = beta of the SNP-exposure association, eaf = effect allele frequency of the SNP, se = standard error of the SNP-exposure association,
### N = sample size of the SNP-exposure GWAS

for(i in 1:length(dat_list)){
  dat <- dat_list[[i]] #harmonization
  eaf <- dat$eaf.exposure     
  b <- dat$beta.exposure
  se <- dat$se.exposure
  p <- dat$pval.exposure
  snp<- dat$SNP
  N <- (dat$samplesize.exposure)
  samplesize <- dat$samplesize.exposure
  Ncontrols <- dat$ncontrol_col
  Ncases <- dat$ncase_col
  k<-1
  R2 = (2*(b^2)*eaf*(1-eaf))/(2*(b^2)*eaf*(1-eaf))+((se^2)*(2*N)*eaf*(1-eaf))
  dat$R2 <- paste0(R2)
  
  #Calculate individual F-stats for this analysis (where F = R2*(n-1-k)/((1-r2)*k), and k is the number of SNPs)
  
  F<-unlist(lapply(1:length(R2),FUN=function(x) R2[x]*(N-1-k)/((1-R2[x])*k)))
  F_stat <- as.character(F)
  dat$F_stat <- F_stat
  rm(F)
  
  #individual power calculation
  
  sig<-0.05 #alpha
  ratio<-Ncontrols/Ncases
  hypothesised_b <- 0.1823216 #ln(OR) of 1.2 in cancer risk 
  power<-pnorm(sqrt(N*R2*(ratio/(1+ratio))*(1/(1+ratio)))*hypothesised_b-qnorm(1-sig/2))
  dat$power <- paste0(power)
  name <- paste0(dat_list[[i]]$exposure[1])
  F_list[[paste(i, "_", dat_list[[i]]$exposure[1], sep = "")]] <- assign(name,dat)
  filename <- paste0("./analysis/", name, "_","harmonised_data",".txt")
  write.table(dat, file = filename, col.names = T, row.names = F, quote = FALSE, sep = "\t")
  rm(name)
}

F_statistics <- rbindlist(F_list)
F_statistics_new <- cbind(F_statistics$exposure, F_statistics$R2, F_statistics$F_stat, F_statistics$power)
names <- c("Exposure", "R2", "F_statistic", "Power")
colnames(F_statistics_new) <- names

write.table(F_statistics, "./analysis/F_R2_Power_MT_CIMBA.txt", col.names = T, row.names = F, quote = FALSE, sep = '\t')

#########################################################################################################
######## PROJECT: Running MR analysis between gut microbial trait abundances/presence and absence and Breast Cancer (from BCAC) with FGFP data
######## Script: Grace Edmunds, adapted from Kaitlin Wade Performing MR analysis with FGFP abundances as exposure 
######## Date: 19/10/20

## clearworkspace
rm(list = ls())

### LOAD PACKAGES 
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

### Prepare environment

## get working directiory
getwd()

## if incorrect - set working directory 
#setwd(<mydirectory>)

#########################################
## PREPARING GUT MICROBIOME GWAS FILES ##
#########################################
#NB: The microbiome traits specified are those that are presented in Table S1 of the paper

#Assign names to run analyses with
bugs <- c("C_Gammaproteobacteria_RNT","F_Sutterellaceae_HB","G_Bifidobacterium_RNT","G_Butyricicoccus_RNT","G_Dialister_HB","G_Parabacteroides_RNT",
          "G_Ruminococcus_HB","G_unclassified_F_Erysipelotrichaceae_HB","G_unclassified_F_Porphyromonadaceae_RNT","G_unclassified_O_Bacteroidales_HB","G_unclassified_P_Firmicutes_HB",
          "G_unclassified_P_Firmicutes_RNT","G_Veillonella_HB","G_Coprococcus_HB")

# Making lists to input results within
results_list <- list()

# Pull out all SNPs that are significant across all exposures
for(i in bugs){
  exposure <- paste("./data/Bugs/",i,".txt",sep="")
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


################################
## READ In BCAC GWAS FILES ##
################################
#NB: GWAS of breast cancer executed by zhang et al. 2020 using data from BCAC consortium 

# I pulled the bug snps from the overall bcac data first using a separate script
cancer <- as.data.frame(fread("./data/overall_bcac_bugsnps.txt"))
#preview it
head(cancer)
colnames(cancer)
#select only the columns we want
cancer_restrict <- cancer [,c(18,9,10,11,12,13,17,19)]
head(cancer)
#change column names 
# colnames(cancer) <- c("rsid", "effect_allele.outcome", "ref_allele.outcome", "beta.outcome", "SE.outcome", "pval.outcome", "samplesize.outcome", "EAF")
dim(cancer_restrict)

cancer_restrict$outcome <- "Overall Breast Cancer"
# write.table(cancer_restrict, "./analysis/201007_outcome_data_overall.txt", row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")

##There is one missing snp, look for a proxy then pull the proxy out of bcac data by returning to the separate script

# missing <- filter(all_exposures, !all_exposures$SNP %in% cancer_restrict$rsid)
##The missing snp is for G_Coprococcus_HB
##LDLink output
# snp_proxy <- read.table("./data/LDlink_coprococcus.txt", header = TRUE, sep="\t")
# dim(snp_proxy)
# snp_best_proxy <- snp_proxy[snp_proxy$R2>=0.90,]
# dim(snp_best_proxy)
# snp_proxies <- as.data.frame(snp_best_proxy$RS_Number)
# colnames(snp_proxies)[1] <- "SNP"
# extended_snplist <- rbind(snp_proxies, all_SNPs)
# extended_snplist <- unique(extended_snplist)
# #remove whitespace 
# extended_snplist2 <- extended_snplist
# extended_snplist <- as.data.frame(gsub(" ", "", extended_snplist$SNP))
# colnames(extended_snplist)[1] <- "SNP"
# extended_snplist <- filter(extended_snplist, !extended_snplist$SNP == ".")

# write.table(extended_snplist, "./data/bugs_and_proxies_for_Overallbc.txt", sep = '\t')

#now go back and pull this list out of the bcac data before continuing
#this proxy is not in the BCAC overall data so we cannot get results for coporococcus

rm(list=ls())

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

# Run the MR analyses
for(i in bugs){
  exposure <- paste("./data/Bugs/",i,".txt",sep="")
  exposure_dat <- read_exposure_data(exposure,sep="\t", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "P_value", eaf_col = "eaf", effect_allele_col = "allele_B", other_allele_col = "allele_A", samplesize_col = "n")
  exposure_dat$exposure <- i
  outcome_dat <- read_outcome_data("./analysis/201007_outcome_data_overall.txt", sep="\t", snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", pval_col = "pval.outcome", eaf_col = "EAF", effect_allele_col = "effect.allele.outcome", other_allele_col = "ref.allele.outcome", phenotype_col = "outcome")
  outcome_dat$samplesize_col <-   247173
  outcome_dat$ncase_col <- 133384 
  outcome_dat$ncontrol_col <- 113789
  dat <- harmonise_data(exposure_dat, outcome_dat, action =2)
  name <- paste0(dat$exposure, dat$outcome)
  filepath = paste0("./analysis/", name, "harmonised_data.txt")
  write.table(dat, filepath, col.names = T, row.names = F, quote = F, sep = '\t')
  mr_results <- mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
  
  # Outputting each MR result into results_list
  a <- cbind.data.frame(mr_results$exposure, mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)
  name <- paste("res_",i,"_","overall_BC",sep="")
  results_list[[paste(i,"_","overall_BC",sep="")]] <- assign(name,a)
  dat_list[[paste(i,"_","overall_BC",sep="")]] <- assign(name,dat)
  rm(a,name)

}

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

# write.table(all_results, "./analysis/220725_mr_results_BCoverall2020_biome.txt", row.names=F, col.names = T, quote = F, sep = "\t")

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
res2$F_stat<- (res2$R2*(res2$samplesize.exposure-1))/((1-res2$R2)*1)
# write.table(res2, "./analysis/R2values_binary_overall.txt", row.names = F, col.names = T, quote = F, sep = '\t')

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
dat$VG <- 2*(b^2)*eaf*(1-eaf)
dat$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
dat$R2 <- dat$VG/dat$PV
dat$Fstat = dat$R2*(N-1-k)/((1-dat$R2)*k)

#individual power calculation
sig<-0.05 #alpha
ratio<-Ncontrols/Ncases
hypothesised_b <- 0.1823216 #ln(OR) of 1.2 in cancer risk 
power<-pnorm(sqrt(N*dat$R2*(ratio/(1+ratio))*(1/(1+ratio)))*hypothesised_b-qnorm(1-sig/2))
dat$power <- paste0(power)
name <- paste0(dat_list[[i]]$exposure[1])
F_list[[paste(i, "_", dat_list[[i]]$exposure[1], sep = "")]] <- assign(name,dat)
filename <- paste0("./analysis/", name, "_","harmonised_data",".txt")
write.table(dat, file = filename, col.names = T, row.names = F, quote = FALSE, sep = "\t")
rm(name)
}

F_statistics <- rbindlist(F_list)
#write.table(F_statistics, "./analysis/F_R2_Power_MT_OverallBC.txt", col.names = T, row.names = F, quote = FALSE, sep = '\t')

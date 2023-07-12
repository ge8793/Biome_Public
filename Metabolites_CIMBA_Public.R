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
getwd()


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
exposure_dat <- read_exposure_data("./data/220512_Shin_Supp5_Formatted.txt", clump = FALSE, sep="\t", snp_col = "SNP", beta_col = "b", se_col = "se", pval_col = "Pval", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", samplesize_col = "n", phenotype_col = "Biochemical")

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

#read in a lookup table to convert bugs rsids to chr_posn and then pull them from BCAC, this is available on request from the BCAC authors
lookup <- as.data.frame(fread("./data/CIMBA_lookuptable.txt"))

#see what columns are called
colnames(lookup)
lookup[1:10,]

#make the column names into the first row
colnames(lookup) <- lookup[1,]
colnames(lookup) <- c("variant_ID", "SNP")

#pull the biome snps out of the lookup
look_filter <- filter(lookup, lookup$SNP %in% exposure_dat$SNP)

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
write.table(merge, "./data/outcome_CIMBA_2020_metabolites.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t') 


#If you would like to use exactly this outcome data, it is available on github at https://github.com/ge8793/Biome_Public/blob/main/outcome_CIMBA_2020_metabolites.txt
#Read in outcome data
merge <- as.data.frame(fread("./data/outcome_CIMBA_2020_metabolites.txt"))

#add sample size as described in 2020 paper methods
merge$ncase <- as.numeric("16592") 
merge$ncontrol <- as.numeric("103901") 
merge$samplesize.outcome <- as.numeric("120493") 

names(merge)

outcome_data <- TwoSampleMR::format_data(merge, type = "outcome", snps = exposure_dat$SNP,
                                         snp_col = "SNP", beta_col = "beta", se_col = "se", 
                                         effect_allele_col = "effect.allele.outcome", other_allele_col = "ref.allele.outcome", 
                                         ncase_col = "ncase", ncontrol_col = "ncontrol", eaf_col = "eaf.outcome", samplesize_col = "samplesize.outcome", pval_col = "pval")


#check unique SNPS in exposure data
SNPS <- as.data.frame(unique(exposure_dat$SNP)) #218 SNPs

#check unique SNPs in outcome data
SNPS <- as.data.frame(unique(outcome_data$SNP)) #214 snps

#harmonise
dat <- harmonise_data(exposure_dat, outcome_data, action =2) #299 obs
check <- unique(dat$SNP) #214 SNPs 
check <- unique(dat$exposure) #217 exposures 

#some will be dropped again now because MR keep = false

mr_results <- mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_results$outcome <- "Triple Negative Breast Cancer"
mr_results$OR <- exp(mr_results$b)
mr_results$LCL <- exp((mr_results$b)-1.96*(mr_results$se))
mr_results$UCL <- exp((mr_results$b)+1.96*(mr_results$se))

check <- unique(mr_results$exposure) #214 unique metabs 

write.table(mr_results, "./analysis/Metabolites/CIMBA_Results/230707_MR_Results_Metabolites_CIMBA.txt", col.names = T, row.names = F, quote = F, sep = '\t')

#get egger p value
het <- mr_pleiotropy_test(dat)

#Get the Q stats
Q <- mr_heterogeneity(dat)
Q_baso <- as.data.frame(Q)
Q$outcome <- paste0("Meta-analysed Triple Negative Breast Cancer")

write.table(Q_baso, "./analysis/Metabolites/Metabolite_Q_statistics_CIMBA.txt", row.names = T, col.names = T, sep = '\t')

#add F stats and R2 if needed (as the exposure data is the same as the other analyses you only need to do this once)

write.table(dat, "./analysis/Metabolites/CIMBA_harmonised_data.txt", col.names = T, row.names = F, quote = F, sep = '\t')

names(exposure_dat)
eaf <- as.numeric(exposure_dat$eaf.exposure)     
b <- as.numeric(exposure_dat$beta.exposure)
se <- as.numeric(exposure_dat$se.exposure)
p <- as.numeric(exposure_dat$pval.exposure)
snp<- exposure_dat$SNP
N <- exposure_dat$samplesize.exposure #number of individuals in shin
samplesize <- exposure_dat$samplesize.exposure
k<-1
exposure_dat$VG <- 2*(b^2)*eaf*(1-eaf)
exposure_dat$VG <- as.numeric(exposure_dat$VG)
exposure_dat$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
exposure_dat$R2 <- exposure_dat$VG/exposure_dat$PV
exposure_dat$Fstat = exposure_dat$R2*(N-1-k)/((1-exposure_dat$R2)*k)
colnames(exposure_dat)
exposure_supp <- exposure_dat[,c(1,2,3,4,5,6,7,8,9,16,17)]
write.table(exposure_supp, "./manuscript/July_2023/Metab_Exposure_Dat_Supp_S3_CIMBA.txt", col.names = T, row.names = F, quote = F, sep = '\t')


#Undertake LOO analysis
loo <- mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)
write.table(loo, "./analysis/Metabolites/CIMBA_Results/CIMBA_LOO.txt", col.names = T, row.names = F, quote = F, sep = '\t')

#Make and export a list of LOO plots
loo_plot <- mr_leaveoneout_plot(loo)
class(loo_plot)
library(RColorBrewer)

pdf(file = "./analysis/Metabolites/CIMBA_Results/LOO_TNBC.pdf")
for (i in 1:length(loo_plot)) {
  print(loo_plot[[i]]+scale_color_manual(values =rep(brewer.pal(5,"Set1"), times=4)))
}
dev.off()



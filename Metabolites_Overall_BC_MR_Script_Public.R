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
## PREPARING EXPOSURE DATA AS METABOLITE GWAS FILES ##
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

#write.table(Shin, "./data/Shin_Supp5_Formatted.txt", sep = '\t', row.names = F, col.names = T, quote = F)

#If you have already done the above, you can read in the exposure data 
exposure_dat <- read_exposure_data("./data/220512_Shin_Supp5_Formatted.txt", clump = FALSE, sep="\t", snp_col = "SNP", beta_col = "b", se_col = "se", pval_col = "Pval", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele", samplesize_col = "n", phenotype_col = "Biochemical")

#or format the Shin data without saving and reading back in
exposure_dat <- TwoSampleMR::format_data(dat = Shin[,-1], type = "exposure",
                                                 phenotype_col = "Biochemical", snp_col = "SNP",
                                                 beta_col = "b", se_col = "se", eaf_col = "eaf",
                                                 effect_allele_col = "effect_allele", 
                                                 other_allele_col = "other_allele", pval_col = "Pval", 
                                                 samplesize_col = "n")

exposure_dat <- exposure_dat[c(2:300),]

#########################################
## PREPARING BCAC DATA FILES AS OUTCOME DATA ##
#########################################

#We downloaded the full BCAC outcome data from https://bcac.ccge.medschl.cam.ac.uk/bcacdata/icogs-complete-summary-results/
#In this analysis we are looking at Overall BC so use the Overall dataset saved as 'icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt'

#The data is large so preview before you read in and then you know you only need columns "SNP.iCOGs", "Effect.Meta", "Baseline.Meta", "Beta.meta", "sdE.meta"
BCAC_preview <- as.data.frame(fread("./data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt", nrows = 10))
names(BCAC_preview)

BCAC <- as.data.frame(fread("./data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt", 
                            select = c("SNP.iCOGs", "Effect.Meta", "Baseline.Meta", "Beta.meta", "sdE.meta", "Freq.Gwas", "EAFcases.iCOGs", "EAFcontrols.iCOGs", "EAFcontrols.Onco", "EAFcases.Onco")))

#calculate the eaf to use 
#weighted average of the contributing studies to the meta-analysis
BCAC$eaf <- (((BCAC$Freq.Gwas * 32498) + (BCAC$EAFcases.iCOGs * 38349) + (BCAC$EAFcontrols.iCOGs * 37818) + (BCAC$EAFcases.Onco * 80125) + (BCAC$EAFcontrols.Onco * 58383)) / 247173)


#see what columns are called
colnames(BCAC)
BCAC[1:10,]

#trim down
BCACtrim <- BCAC[,c(1,2,3,4,5,11)]
names(BCACtrim)
colnames(BCACtrim) <- c("snp.icogs", "effect.allele.outcome", "ref.allele.outcome", "beta.outcome", "se.outcome", "eaf.outcome")

#add sample size as described in 2020 paper methods supplementary table 4
BCACtrim$ncase <- as.numeric("133384") 
BCACtrim$ncontrol <- as.numeric("113789") 
BCACtrim$samplesize.outcome <- as.numeric("247173") 

head(BCACtrim)

#split the snp column to get just the rsid out of it 
BCACtrim[,10:13] <- str_split_fixed(BCACtrim$snp.icogs, ":", 4)
#preview the first ten rows
BCACtrim [1:10,]
#rename the rsid column
colnames(BCACtrim)[10] <- "rsid"
#preview the first ten rows
BCACtrim [1:10,]
#delete the last 3 columns
BCACtrim <- BCACtrim[,1:10]
#preview the first ten rows
BCACtrim [1:10,]

#NB there is no p-value column, infer from MR base, do not use the pvalues in the BCAC data, they are not really the GWAS pvalue

#output this new dataframe as a csv file
write.csv(BCACtrim, "./data/BCAC_overall_2020_trimmed.csv", sep = '/t', row.names = FALSE, col.names = TRUE, quote = FALSE) 

## Format overall breast cancer outcome data 

names(BCACtrim)

outcome_data <- TwoSampleMR::format_data(BCACtrim, type = "outcome", snps = exposure_dat$SNP,
                                    snp_col = "rsid", beta_col = "beta.outcome", se_col = "se.outcome", 
                                    effect_allele_col = "effect.allele.outcome", other_allele_col = "ref.allele.outcome", 
                                    ncase_col = "ncase", ncontrol_col = "ncontrol", eaf_col = "eaf.outcome")

outcome_data[1:10,]

write.table(outcome_data, "./data/Outcome_overall_BCAC_metabolites.txt", row.names = F, col.names = T, quote = F, sep = '\t')

#If you want to use exactly this outcome data you can download it from Github available at https://github.com/ge8793/Biome_Public/blob/main/outcome_overall_BCAC_metabolites.txt
outcome_data <- as.data.frame(fread("./data/Outcome_overall_BCAC_metabolites.txt"))

## N.B. MR-Base infers p-values from other fields because p-values are not provided

###### Prep for MR

#Optional check load Ryan outcome data and run with this instead
#load("./data/overall_outcome.Robj")
#outcome_data <- overall

#check unique SNPS in exposure data
SNPS <- as.data.frame(unique(exposure_dat$SNP)) #218 SNPs
exposures <- as.data.frame(unique(exposure_dat$exposure)) 
#check unique SNPs in outcome data
SNPS <- as.data.frame(unique(outcome_data$SNP)) #204 SNPs

#Make a list of the SNPs that are missing and manually look them up in LDlink
proxy_list <- filter(exposure_dat, !exposure_dat$SNP %in% outcome_data$SNP) #19 snps
missing <- unique(proxy_list$SNP) #14 exposures 

#harmonise
dat <- harmonise_data(exposure_dat, outcome_data, action =2)
check <- unique(dat$SNP) #204 SNPs
check <- unique(dat$exposure) #209 exposures

#some will be dropped again now because MR keep = false

mr_results <- mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_results$outcome <- "OverallBC"
mr_results$OR <- exp(mr_results$b)
mr_results$LCL <- exp((mr_results$b)-1.96*(mr_results$se))
mr_results$UCL <- exp((mr_results$b)+1.96*(mr_results$se))

write.table(mr_results, "./analysis/Metabolites/Overall_Results/230705_MR_Results_Metabolites_OverallBC.txt", col.names = T, row.names = F, quote = F, sep = '\t')

#mr_results <- as.data.frame(fread("./analysis/Metabolites/Overall_Results/MR_Results_Metabolites_OverallBC.txt"))
Res <- mr_results %>% filter(mr_results$method %in% c("Wald ratio", "Inverse variance weighted"))
count <- unique(Res$id.exposure)

#get egger p value
het <- mr_pleiotropy_test(dat)

#Get the Q stats
Q <- mr_heterogeneity(dat)
Q_baso <- as.data.frame(Q)
Q$outcome <- paste0("Overall Breast Cancer")

write.table(Q_baso, "./analysis/Metabolites/Overall_Results/Metabolite_Q_statistics_Overall.txt", row.names = T, col.names = T, sep = '\t')

#add F stats and R2 to harmonised data
names(dat)
eaf <- as.numeric(dat$eaf.exposure)     
b <- as.numeric(dat$beta.exposure)
se <- as.numeric(dat$se.exposure)
p <- as.numeric(dat$pval.exposure)
snp<- dat$SNP
N <- dat$samplesize.exposure #number of individuals in shin
samplesize <- dat$samplesize.exposure
k<-1
dat$VG <- 2*(b^2)*eaf*(1-eaf)
dat$VG <- as.numeric(dat$VG)
dat$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
dat$R2 <- dat$VG/dat$PV
dat$Fstat = dat$R2*(N-1-k)/((1-dat$R2)*k)

write.table(dat, "./analysis/Metabolites/Overall_Results/Overall_harmonised_data.txt", col.names = T, row.names = F, quote = F, sep = '\t')

#add F stats and R2 to exposure data to make supplementary tables 
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
write.table(exposure_supp, "./manuscript/July_2023/Metab_Exposure_Dat_Supp_S3.txt", col.names = T, row.names = F, quote = F, sep = '\t')


#Undertake LOO analysis
loo <- mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)
write.table(loo, "./analysis/Metabolites/Overall_Results/Overall_LOO.txt", col.names = T, row.names = F, quote = F, sep = '\t')

#Make and export a list of LOO plots
loo_plot <- mr_leaveoneout_plot(loo)
class(loo_plot)
library(RColorBrewer)

pdf(file = "./analysis/Metabolites/Overall_Results/LOO_Overall.pdf")
for (i in 1:length(loo_plot)) {
  print(loo_plot[[i]]+scale_color_manual(values =rep(brewer.pal(5,"Set1"), times=4)))
}
dev.off()

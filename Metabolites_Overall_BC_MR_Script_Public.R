#########################################################################################################
######## PROJECT: Running MR analysis between metabolites and Breast Cancer (from BCAC) with Shin et al. metabolite data
######## Script: Grace Edmunds
######## Date: 27/01/22


###### PART 1 - Pull metabolite SNPS out as exposure data

### LOAD PACKAGES
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
library(stringr)
ga_auth(email="graceledmunds@gmail.com")
library(MRInstruments)
library(fuzzyjoin)

#clearws
rm(list=ls())

#Read in the list of FGFP metabolites

mets <- as.data.frame(fread("./data/FGFP_2020_05_21_QCd_feature_data_new_metabolite_list.txt"))
mets$name <- mets$BIOCHEMICAL

# Extract exposure instruments from Shin et al. in IEU GWAS catalogue ====
data(gwas_catalog)
head(gwas_catalog)

mets_gwas <-
  subset(gwas_catalog,
         grepl("metabolite", Phenotype))

table(mets_gwas$Author)

shin_gwas <-
  subset(gwas_catalog,
         grepl("Shin", Author))

shin_gwas$name <- shin_gwas$Phenotype_info

#fuzzyjoin to find the metabolites we need
joined <- stringdist_inner_join(shin_gwas, mets, by = "name")
names(joined)
#run a quick check
check <- as.data.frame(cbind(joined$name.x, joined$name.y))
head(check)

#create exposure data
exposure_data <- filter(shin_gwas, Phenotype_info %in% check$V1)

#Make suplicate exposure data with COMP _IDs attached as will need this info later
names(joined)
comp_IDs <- joined[,c(8,12,29,32)]
names(comp_IDs)
exposure_data2 <-stringdist_inner_join(exposure_data, comp_IDs, by = "Phenotype_info")
check2 <- as.data.frame(cbind(exposure_data2$Phenotype_info.x, exposure_data2$Phenotype_info.y))
names(exposure_data2)
exposure_data2 <-exposure_data2[,c(1:26, 29, 30)]
names(exposure_data2)
names(exposure_data2)[12] <- "SNP"
exposure_data2 <-unique(exposure_data2)
# write.table(exposure_data, "./data/draft002/COMP_IDs_FGFP_metabolites_from_shin.txt", quote = F, row.names = F, col.names = T, sep = '\t')

names(exposure_data)
exposure_data <- format_data(exposure_data, phenotype_col = "Phenotype_info")
met_list <- as.data.frame(table(exposure_data$exposure))
length(met_list$Var1)

#write.table(exposure_data, "./data/exp_dat_FGFP_metabolites_from_shin.txt", quote = F, row.names = F, col.names = T, sep = '\t')
exposure_data <- as.data.frame(fread("./data/exp_dat_FGFP_metabolites_from_shin.txt"))

###### PART 2 - PULL THE SNPS OUT OF BCAC AS OUTCOME DATA ############################################

### Overall BC

#read in metabolite exposure data
metabolites <- as.data.frame(fread("./data/exp_dat_FGFP_metabolites_from_shin.txt"))
head(metabolites)
SNPs <- as.data.frame(table(metabolites$SNP)) #This tells us there are 72 individual SNPs in the data 

#read in BCAC Overall outcome data 
check <- as.data.frame(fread("./data/overall_2020.csv"))
check[1:10,]

#filter on the exposure SNPs
bcac_filter <- filter(check, check$rsid %in% metabolites$SNP) #69/72 SNPs were available in BCAC Overall BC

#Make a list of the SNPs that are missing and manually look them up in LDlink, add any proxies into the data 
proxy_list <- filter(metabolites, !metabolites$SNP %in% bcac_filter$rsid)

#E.g. for rs1498694
snp_proxy <- read.table("./data/Metabolite_Proxies/LDLINK_rs1498694.txt", header = TRUE, sep="\t")
dim(snp_proxy)
snp_best_proxy <- snp_proxy[snp_proxy$R2>=0.90,]
dim(snp_best_proxy)
snp_proxies <- as.vector(snp_best_proxy$RS_Number)
test <- check[check$rsid %in% snp_proxies,] #1 proxies available rs1722383
test$original <- "rs1498694"
test$proxy <- "1"
rs1498694_dat <- test

#Add rs1498694 proxy to the outcome data
bcac_filter$original <- bcac_filter$rsid
bcac_filter$proxy <- "0"
bcac_extended <- rbind.data.frame(bcac_filter, rs1498694_dat)

## save data
# write.table(bcac_extended, "./analysis/Metabolites/outcome_data_overall_metabolites_with_proxies.txt", row.names = F, col.names = T, quote = F, sep = '\t')

#Remove all of this this now as it is big
rm(list=ls())

###### PART 3 - RUN MR

### Overall BC

## Read in metabolite exposure data and continue formatting

# metabolites <- as.data.frame(fread("./data/exp_dat_FGFP_metabolites_from_shin.txt"))
# head(metabolites)

# exposure_dat <- metabolites
# head(exposure_dat)
# 
# exposure_dat$exposurenew <- str_replace_all(exposure_dat$exposure, "\\(unit increase\\)", "")
# exposure_dat$exposurenew <- str_replace_all(exposure_dat$exposurenew, "\\(unit decrease\\)", "")
# exposure_dat$exposurenew <- str_replace(exposure_dat$exposurenew, "[)]", "")
# exposure_dat$exposurenew <- str_replace(exposure_dat$exposurenew, "[)]", "")
# exposure_dat$exposurenew <- str_replace(exposure_dat$exposurenew, "[(]", "")
# exposure_dat$exposurenew <- str_replace(exposure_dat$exposurenew, "[(]", "")
# exposure_dat$exposure <- exposure_dat$exposurenew
# head(exposure_dat)

# write.table(exposure_dat, "./data/exp_dat_formatted_shin.txt", row.names = F, col.names = T, quote = F, sep = '\t')

exposure_dat <- read_exposure_data("./data/exp_dat_formatted_shin.txt", clump = T, sep="\t", snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure", pval_col = "pval.exposure", eaf_col = "eaf.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", phenotype_col = "exposure")

## Read in BCAC Overall outcome data and continue formatting
cancer <- as.data.frame(fread("./analysis/Metabolites/outcome_data_overall_metabolites_with_proxies.txt"))

#make a copy of the data in case anything goes wrong
#cancer2 <- cancer
#preview it
#head(cancer)
#colnames(cancer)
#select only the columns we want
#cancer_restrict <- cancer [,c(18,9,10,11,12,13,17,19)]
#head(cancer_restrict)
#change column names 
#colnames(cancer_restrict) <- c("rsid", "effect_allele.outcome", "ref_allele.outcome", "beta.outcome", "SE.outcome", "pval.outcome", "samplesize.outcome", "EAF")
#dim(cancer_restrict)
# write.table(cancer_restrict, "./analysis/Metabolites/Overall_outcome_data_fotmatted.txt", col.names = T, row.names = F, quote = F,  sep = '\t')

outcome_dat <- read_outcome_data("./analysis/Metabolites/Overall_outcome_data_fotmatted.txt", sep="\t", snp_col = "rsid", beta_col = "beta.outcome", se_col = "SE.outcome", pval_col = "pval.outcome", eaf_col = "EAF", effect_allele_col = "effect_allele.outcome", other_allele_col = "ref_allele.outcome", phenotype_col = "Overall_Breast_Cancer")
outcome_dat$samplesize_col <-   247173
outcome_dat$ncase_col <- 133384 
outcome_dat$ncontrol_col <- 113789
outcome_dat$Phenotype <- "Overall_BC"

dat <- harmonise_data(exposure_dat, outcome_dat, action =2)

names(dat)

#Get the MR Egger intercept for any with >2 SNPs 
het <- mr_pleiotropy_test(dat)

mr_results <- mr(dat, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_results$outcome <- "OverallBC"
mr_results$OR <- exp(mr_results$b)
mr_results$LCL <- exp((mr_results$b)-1.96*(mr_results$se))
mr_results$UCL <- exp((mr_results$b)+1.96*(mr_results$se))
#write.table(mr_results, "./analysis/Metabolites/overall_results.txt", col.names = T, row.names = F, quote = F, sep = '\t')

#Get the Q stat for any with >2 SNPS
Q_list <- list ()
Q <- mr_heterogeneity(dat)
Q_baso <- as.data.frame(Q)
Q$outcome <- paste0("Overall Breast Cancer")
#write.table(Q_baso, "./analysis/Metabolites/Metabolite_Q_statistics_Overall.txt", row.names = T, col.names = T, sep = '\t')

#Get the R2 and F stat for each SNP

eaf <- as.numeric(dat$eaf.exposure)     
b <- as.numeric(dat$beta.exposure)
se <- as.numeric(dat$se.exposure)
p <- as.numeric(dat$pval.exposure)
snp<- dat$SNP
N <- 7824 #number of individuals in shin
samplesize <- dat$samplesize.exposure
Ncontrols <- dat$ncontrol_col
Ncases <- dat$ncase_col

dat$VG <- 2*(b^2)*eaf*(1-eaf)
dat$VG <- as.numeric(dat$VG)
dat$PV <- 2*(b^2)*eaf*(1-eaf) + ((se^2)*2*N*eaf*(1-eaf))
dat$R2 <- dat$VG/dat$PV

##The F-stat will be different for single and multi-SNP instruments 

#single
single <- filter(mr_results, mr_results$method == "Wald ratio")
dat_single <- filter(dat, dat$exposure %in% single$exposure) 
k <- 1 #nsnp
dat_single$Fstat = dat_single$R2*(N-1-k)/((1-dat_single$R2)*k)

#multi
multi <- filter(mr_results, !mr_results$method == "Wald ratio")
names(multi)
multi_trim <- multi[,c(4,6)]
dat_multi <- filter(dat, dat$exposure %in% multi$exposure) 
dat_multi <- merge(dat_multi, multi_trim, by= "exposure")
k <- dat_multi$nsnp #nsnp
dat_multi$Fstat = dat_multi$R2*(N-1-k)/((1-dat_multi$R2)*k)

#write.table(dat_single, "./analysis/Metabolites/Overall_harmonised_data_singlesnps.txt", col.names = T, row.names = F, quote = F, sep = '\t')
#write.table(dat_multi, "./analysis/Metabolites/Overall_harmonised_data_multisnps.txt", col.names = T, row.names = F, quote = F, sep = '\t')


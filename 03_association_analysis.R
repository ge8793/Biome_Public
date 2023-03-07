###########################################
## Asociation analysis for Kaitlin
## - FGFP Microbiome on metabolome
## - By: David Hughes
## - Date: March 18th 2022
###########################################


###########################################
## prep
###########################################
library(tidyverse)

data_dir = "/group/njt-group/fgfp_4_KW/data/"
source("rntransform.R")
source("lmfit.R")
source("normW.R")

###########################################
## Read in and format the data
###########################################
## covariables
f = paste0(data_dir, "my_covariable_data.txt")
cov_data = read.table(f, header = TRUE, sep = "\t", as.is = TRUE)

###########################################
## Edit Variables
###########################################
## turn sex and is_vegetarian into a factor
w = c(2,9); for(i in w){ cov_data[,i] = as.factor(cov_data[,i]) }
## Make smoker never|ever|current into a numeric
cov_data$smoker_NEC[ cov_data$smoker_NEC == "never"] = 0
cov_data$smoker_NEC[ cov_data$smoker_NEC == "ever"] = 1
cov_data$smoker_NEC[ cov_data$smoker_NEC == "current"] = 2
cov_data$smoker_NEC = as.numeric(cov_data$smoker_NEC)
## Make meat consumption numeric
cov_data$meat_consumption[ cov_data$meat_consumption == "Not applicable"] = 0
cov_data$meat_consumption[ cov_data$meat_consumption == "Between 5 and 10"] = 6
cov_data$meat_consumption[ cov_data$meat_consumption == "More than 10"] = 7
cov_data$meat_consumption = as.numeric(cov_data$meat_consumption)
## Make meat frequency numeric
cov_data$meat_frequency[ cov_data$meat_frequency == "I hardly ever or never eat/drink these products"] = 0
cov_data$meat_frequency[ cov_data$meat_frequency == "Not this week"] = 1
cov_data$meat_frequency[ cov_data$meat_frequency == "One time"] = 2
cov_data$meat_frequency[ cov_data$meat_frequency == "A few times"] = 3
cov_data$meat_frequency[ cov_data$meat_frequency == "Almost every day"] = 4
cov_data$meat_frequency[ cov_data$meat_frequency == "Every day"] = 5
cov_data$meat_frequency[ cov_data$meat_frequency == "Multiple times a day"] = 6
cov_data$meat_frequency = as.numeric(cov_data$meat_frequency)


###########################################
## microbial traits data
###########################################
f = paste0(data_dir, "my_microbial_traits_data.txt")
mt_data = read.table(f, header = TRUE, sep = "\t", as.is = TRUE)

## MAKE PA MTs a factor
w = grep("PA", colnames(mt_data) )
for(i in w){
  mt_data[,i] = as.factor(mt_data[,i])
}
## Set Enterotype as a factor
mt_data$EnterotypeClass = as.factor(mt_data$EnterotypeClass)

## A vector of microbial trait IDs
micro_trait_ids = colnames(mt_data)[-c(1:10)]

## MT phylogeny data
f = paste0(data_dir, "microbiome_phylogeny.txt")
phylogeny_data = read.table(f, header = TRUE, sep = "\t", as.is = TRUE)

## metabolite data
f = paste0(data_dir, "metabolite_data.txt")
metabolite_data = read.table(f, header = TRUE, sep = "\t", as.is = TRUE)
## filter those with less than 1000 observations
not_na_count = apply(metabolite_data, 2, function(x){ sum(!is.na(x)) })
w = which(not_na_count < 1000)
metabolite_data = metabolite_data[,-w]

## A vector of metabolite IDs
metabo_ids = colnames(metabolite_data)

## metabolite feature data
f = paste0(data_dir, "metabolite_feature_data.txt")
f = metabolite_feature_file
feature_data = read.table(f, header = TRUE, 
                          sep = "\t", 
                          as.is = TRUE,
                          quote="")


###########################################
## Define mydata
###########################################
mydata = cbind(cov_data, mt_data, metabolite_data )


###########################################
## Run the association analysis
###########################################
SumData = c()

for(mt in micro_trait_ids[ 1:length(micro_trait_ids) ] ){
  cat(paste0("now running ", mt, "\n"))
  for(metabolite in metabo_ids){
    o = lmfit( wdata = mydata,
               outcome = metabolite,
               exposure = mt,
               covariates = c("sex","age","bmi",
                              "smoker_NEC","drinking_week_avg",
                              "is_vegetarian","monthly_hh_income"),
               weights = NA,
               rnt_outcome = TRUE,
               typeIIanova = TRUE)
    ###
    SumData = rbind(SumData, o)
  }
}
  
###
SumData = as.data.frame(SumData)
for(i in 3:ncol(SumData)){ SumData[,i] = as.numeric(SumData[,i]) }


###########################################
## Add Metabolite Annotation to the results
###########################################

m = match(SumData$outcome, feature_data$feature_names)
SumData = cbind(SumData, feature_data[ m , c("BIOCHEMICAL", "CHEMICAL_ID", "HMDB", "SUPER_PATHWAY", "SUB_PATHWAY", "PUBCHEM", "KEGG" )] )

###########################################
## Write data to file
###########################################
f = paste0("/group/njt-group/fgfp_4_KW/results/Association_SumStats_v0.2.txt")

write.table(SumData, file = f, 
            row.names = FALSE, col.names = TRUE, 
            sep = "\t", quote = FALSE)










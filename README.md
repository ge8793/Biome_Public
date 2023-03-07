# FGFP: microbial trait metabolite association analysis

### - breast cancer study  for Kaitlin Wade, Tim Robinson, and Grace Edmunds 


### - request from Kaitlin

Hey moose! Would it be possible to run those regression analyses between the microbiota and metabolites for the breast cancer paper? Well only need three microbiota and I can pick those out if you want but if its easy to just do a loop of everything then we can pick out those results we need. I think we decided on the following variables:

Outcomes: - All metabs

Exposures: - All microbiota measures (same ones we used for GWAS

Confounders: 
	
	age (age)
	sex (gender)
	bmi (BMI)
	smoking (have_smoked_foafull_year)
	alcohol (alcohol_average_consumption_last_week)
	marker of vegetarian diet 
		e.g. could make a categorical variable from the "is_vegan" and "is_vegetarian" variables with those who haven't specified either we can assume are meat eaters
	marker of socioeconomic position 
		do we have anything like this? education? income?
		
## Model

	fit = lm(metabolite ~ sex + age + bmi + smoker_NEC + drinking_week_avg + is_vegetarian + monthly_hh_income + micobiome_trait)
	
	- the outcomes were rank normal transformed

## NOTE:
	- There were some obvious errors in the anthropomorphic triats in the newest release of FGFP. 
	- As such, on March 17th 2022 David H. prepared a new release of the FGFP clinical data file and re-ran this analysis. Result file is now "../results/Association_SumStats_v0.2.txt"

## variables written to file

	1. exposure: the names of the exposure (microbial trait)
	2. outcome: the name of the outcome (metabolite)
	3. n: the sample size of this analysis
	4. W: the Shapiro-Wilk (normality) W-statistic for the model residuals
	5. Breusch_Pagan_P: The P-value for Breusch-Pagan heteroskedasticity test
	6. beta: the effect estimate of exposure on outcome (normalized standard deviation units per 1 unit increase (integer count) in exposure)
	7. se: standard error of beta
	8. tval: t-statistic for exposure outcome association t-test
	9. P: p-value for exposure outcome association
	10. sw_beta: exposure outcome beta (again) but to accompany the robust sandwich SE estimates
	11. sw_se: robust sandwich estimate of the standard error (SE) in the beta. variance covariance matrix estiamted with White’s estimator
	12. sw_zval: z-statistic for exposure outcome (robust sandwich) association
	13. sw_P: p-value for robust sandwich SE exposure outcome association z-test
	14. model_R2: the total variance in outcome explained by the model
	15. model_adjR2: an adjusted estimate of the total variance in outcome explained by the model. Only useful if number of model predictors change.
	16. etasq_sex: the variance in outcome explained by sex (derived type II ANOVA sums of squares)
	17. etasq_age: the variance in outcome explained by age (derived type II ANOVA sums of squares)
	18. etasq_bmi: the variance in outcome explained by body mass index (derived type II ANOVA sums of squares)
	19. etasq_smoker_NEC: the variance in outcome explained by smoking status [additive Never-Ever-Current] (derived type II ANOVA sums of squares)
	20. etasq_drinking_week_avg: the variance in outcome explained by average weekly alcohol consumption (derived type II ANOVA sums of squares)
	21. etasq_is_vegetarian: the variance in outcome explained by vegetarians (derived type II ANOVA sums of squares)
	22. etasq_monthly_hh_income: the variance in outcome explained by monthly household income (derived type II ANOVA sums of squares)
	23. etasq_exposure: the variance in outcome explained by exposure (microbial trait) (derived type II ANOVA sums of squares)
	24. etasq_Residuals: the unexplained (residual | error) variance in outcome (derived type II ANOVA sums of squares)
	25. sex_P: type II anova p-value for sex (same P as in linear model t-test)
	26. age_P: type II anova p-value for age (same P as in linear model t-test)
	27. bmi_P: type II anova p-value for body mass index (same P as in linear model t-test)
	28. smoker_NEC_P: type II anova p-value for smoking status [additive: Never-Ever-Current] (same P as in linear model t-test)
	29. drinking_week_avg_P: type II anova p-value for average weekly alcohol consumption (same P as in linear model t-test)
	30. is_vegetarian_P: type II anova p-value for being a vegetarian (same P as in linear model t-test)
	31. monthly_hh_income_P: type II anova p-value for monthly household income (same P as in linear model t-test)
	32. exposure_P: type II anova p-value for exposure (same P as in linear model t-test)
	33. BIOCHEMICAL: metabolon provided biocemical ID or metabolite name
	34. CHEMICAL_ID: metabolon provided chemical ID
	35. HMDB: metabolon provided Human Metabolome Database ID
	36. SUPER_PATHWAY: metabolon provided super pathway category
	37. SUB_PATHWAY: metabolon provided sub pathway category
	38. PUBCHEM: metabolon provided PubChem ID
	39. KEGG: metabolon provided KEGG pathway ID

	
## observed associations 

### A tile plot of associations with a P < 0.05/228682

![heatmap](figures/log10_P_tileplot_v0.2.png)


# This is the header file containing functions for the simulation in Senior et al. Bias in meta-analysis of response ratios is reduced by preliminary meta-analysis of the variance. 

# This script contains 4 functions that are used at various points in the proceeding sets of scripts.

# The first function listed is sim_data
# The second function listed is my_lnRR
# The third function listed is my_lnCV
# The fourth function listed is my_median

# Comments and descriptions for each preceed the code.

############################################
############### sim_data ###################
############################################

# The function is designed to simulate multi-level meta-analytic data. The are simulated hierarchically at the among-study and within-study levels.

# The same function however can be used to simulate single-level data (e.g., a conventional random-effects meta-analysis, where there is just one effect per study) by setting k_effect_mu at 1, k_effect_sd at 0, k_study at the desired number of studies, and icc_study = 1. This parameterisation forces all the heterogeneity at the study level, and there to be no within-study effect sizes.

# The function allows the user to simulate heterogeneous (among-study) differences in sample size and among-sample variances. Note that variances and sample-sizes vary at the level of study rather than effect size (or treatment group). That is all groups/samples within the same study will have the same n and SD.

# The arguments to the function are as follows:
# lnRR = the overall mean effect size - assumed to be normally distirbuted (i.e., the RR is log-normal)
# k_study = the number of studies - specifed by the user and in-variant
# k_effect_mu = the mean number of effect sizes per study - the number of effect sizes is double poisson distributed
# k_effect_sd = the sd in number of effect sizes per study - the number of effect sizes is double poisson distributed
# tau2 = the total heterogeneity among effect sizes
# icc_study = the proportion (intraclass correlation) of tau2 coming from the study level effects 
# n_study_mu = the mean sample size studies in the dataset - the sample size is double poisson distributed
# n_study_sd = the sd of sample size of studies in the dataset - the sample size is double poisson distributed
# sd_study_mu is the mean sd of data within effect size for each study - the sds for samples are assumed to be gamma distributed
# sd_study_sd is the sd of sd of data within effect size for each study - the sds for samples are assumed to be gamma distributed
# return_true_effects allows the user to see the simulated effect at the level of each study and effect size in the outputted data

# Returns a dataframe of sample means, sds and ns for the studies, which can then be used to calculate effect sizes.

sim_data<-function(lnRR, k_study, k_effect_mu, k_effect_sd, tau2, icc_study, n_study_mu, n_study_sd, sd_study_mu, sd_study_sd, return_true_effects=F){
	
	# Load the package for the double-poisson distribution
	require(gamlss.dist)
	
	# Start by reparameterizing a few things.
	
	# Make sure the various sds for the sample size and sample-sd distributions are not <= 0, as this is undefined for the gamma and double poisson distributions
	if(sd_study_sd <= 0){
		sd_study_sd<-10^-8
	}
	if(k_effect_sd <= 0){
		k_effect_sd<-10^-8
	}
	if(n_study_sd <= 0){
		n_study_sd<-10^-8
	}
	
	# Calculate the between and within-study heterogeneity from the tau2 and icc study
	study_sigma2<-tau2 * icc_study
	resid_sigma2<-tau2 * (1 - icc_study)
	
	# For the gamma distirbution for the variance the a (shape) and b (scale) are approximated as
	a<-sd_study_mu^2 / sd_study_sd^2
	b<-sd_study_sd^2 / sd_study_mu
	
	# For the double poisson distribution the sigma parameter as follows. I will specify the mean as n_study_mu - 3, then add 3 to avoid 0 sample sizes - same for k_effects. Note the addition of 10^-8 avoids 0 means
	sigma_n<-n_study_sd^2 / (n_study_mu - 3 + 10^-8)
	sigma_k<-k_effect_sd^2 / (k_effect_mu - 1 + 10^-8)
 	
 	# Loop to do the simulations for the k studies
	for(i in 1:k_study){
		
		# Get the study specific-effect for study i
		ES_i<-rnorm(1, lnRR, sqrt(study_sigma2))
		
		# Find study-specific number of effect sizes, sample size and variance for study i
		k_effect_i<-rDPO(1, mu=(k_effect_mu - 1 + 10^-8), sigma=sigma_k) + 1
		n_i<-rDPO(1, mu=(n_study_mu - 3 + 10^-8), sigma=sigma_n) + 3
		sd_i<-rgamma(1, shape=a, scale=b)
		
		# Loop to do the simulations the k effects
		for(j in 1:k_effect_i){
			
			# Get the within-study effects for study i effect size j
			ES_ij<-rnorm(1, ES_i, sqrt(resid_sigma2))	
		
			# Simulate the control group for study i, effect size j
			control_ij<-rnorm(n=n_i, mean=mu_control, sd=sd_i)
			
			# Simulate the treatment group for study i, effect size j
			treatment_ij<-rnorm(n_i, mean(mu_control * exp(ES_ij)), sd=sd_i)
		
			# Now calculate the mean, sd and give n in a nice way for downstream code
			data_ij<-data.frame(Effect = paste(i, j, sep="_"), Study = i, Control.Mean = mean(control_ij), Control.SD = sd(control_ij), Control.n = length(control_ij), Treatment.Mean = mean(treatment_ij), Treatment.SD = sd(treatment_ij), Treatment.n = length(treatment_ij))
			
			# If you want the true simulated effects returned, add those in
			if(return_true_effects){
				data_ij$ES_i<-ES_i
				data_ij$ES_ij<-ES_ij
				data_ij$sd_i<-sd_i
			}
			
			# On the first pass create a data variable from the ijth simulaiton, otherwise bind them all together
			if(i == 1 & j == 1){
				data<-data_ij
			}else{
				data<-rbind(data, data_ij)	
			}	
		}
	}
	
	# Return the data
	return(data)
	
}

############################################
################ my_lnRR ###################
############################################

# Function to calcualte the lnRR based on coefficient of variation (CV) as given in equations 1 and 2 of the manuscript.

# The arguments to the function are as follows:
# Control.Mean = the mean of the control group
# Treatment.Mean = the mean of the treatment group
# Control.CV = the CV of the control group
# Treatment.CV = the CV of the treatment group
# Control.n = the sample size of the control group
# Treatment.n = the sample size of the treatment group

# Note that in the case of effects sizes using the pre-meta-analysis method, the CVs could be estimates, rather than sample specific CVs

# The function returns a vector of two values comprising the estimate of the lnRR and its sampling variance. 

my_lnRR<-function(Control.Mean, Treatment.Mean, Control.CV, Treatment.CV, Control.n, Treatment.n){
	
	# The Effect size
	yi<-log(Treatment.Mean / Control.Mean) + 0.5 * (Treatment.CV^2 / Treatment.n - Control.CV^2 / Control.n)

	# The sampling variance
	vi<-(Treatment.CV^2 / Treatment.n) + (Control.CV^2 / Control.n) + (Treatment.CV^4 / (2 * Treatment.n^2)) + (Control.CV^4 / (2 * Control.n^2))
	
	# Return the data
	return(data.frame(yi=yi, vi=vi))
	
}

############################################
################ my_lnCV ###################
############################################

# Function to calcualte the lnCV and its sampling variance for a single sample based on CV as given in equations 6 and 7 in the manuscript for pre-meta-analysis

# The function takes the following arguments:
# Mean = the mean of the smaple
# SD = the standard deviation of the sample
# n = the sample size of the sample

# The function returns a vector of two values comprising the estimate of the lnCV and its sampling variance. 

my_lnCV<-function(Mean, SD, n){
	
	# The Effect size
	yi<-log(SD / Mean) + (1/(2*(n-1))) + (SD^2/(n * Mean^2))

	# The sampling variance
	vi<-(SD^2/(n * Mean^2)) + (SD^4 / (2 * n^2 * Mean^4)) + (1/(2 * (n - 1)))
 	
	# Return the data
	return(data.frame(yi=yi, vi=vi))
	
}

############################################
################ my_median #################
############################################

# Function to calcualte the median and 99% CI of a set of data. The function uses the following specifcation
# https://www.statology.org/confidence-interval-for-median/

# The function takes just one argument:
# x = the input data

# returns a vector of three values comprsing the lower confidence limit, the median, and the upper confidence limit.

my_median<-function(x){
	
	# Sort the data and get the sample sizes	
	x<-sort(x)
	n<-length(x)
	
	# Find the median
	med<-median(x)
	
	# Find the location of the lower CI value and return it.
	tag<-n*0.5 - 2.58*sqrt(n * 0.5^2)
	lower<-x[ceiling(tag)]
	
	# Find the location of the lower CI value and return it.
	tag<-n*0.5 + 2.58*sqrt(n * 0.5^2)
	upper<-x[ceiling(tag)]
	
	# Return the vector of values
	return(c(lower, med, upper))
	
}

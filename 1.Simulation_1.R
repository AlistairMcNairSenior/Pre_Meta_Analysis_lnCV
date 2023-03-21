
# Simulation to support new analysis looking at step 1 lnCV to improve analysis of lnRR

# This script runs the simulation for the MLMA (i.e. multilevel meta-analysis). It is designed to run as an array job on the HPC and Sydney UNI. There are only 10 reps in here, but I ran on an array 0-999

# Clean up the R Environment 
rm(list=ls())

# Incoming arguments from bash
args<-commandArgs()

# Where are we working
directory<-"/project/RDS-FSC-EvolNutStrats-RW/Miss_Sim"
setwd(directory)

# Load the relevant libraries, and header
library(metafor)
library(plyr)
source("0.Header.R")

# Get the incoming info on the PBS array index
index<-args[6]

# In this simulation I will include not non-independence. The the analysis will be REMA.

###################################################
################## Parameters #####################
###################################################

# Specify the overall effect magnitude - lets assume a lnRR of 0.3 - this is effectively a 35% increase in the mean
lnRR<-0.3

# Specify the number of effect sizes per study 1 with SD 0 for REMA
k_effect_mu<-1
k_effect_sd<-0

# Specify total tau2 of lnRR^2. This sets 0 at 1 SD below the mean effect meaning there will be a set proportion of negative effects. A large degree of heteorgeneity
tau2<-(lnRR * c(0.1, 1))^2

# Specify the % tau2 - 0 and 50%
icc_study<-c(1)

# Note I assume the mean in the control group is 100 - to give +ve means for lnRR
mu_control<-100

# Specify the mean sample size, and its variance - to be simulated from a double poisson distribution
n_study_mu<-c(10, 30)
n_study_sd<-c(sqrt(1.5 * 10), sqrt(1.5 * 30))
# Note we don't want all combos - I'll edit out the low SD with low mean and vice versa

# Specify the mean study SD in outcome, I will set to 15% of mu_control; i.e. the CV is 0.15
sd_study_mu<-mu_control * 0.15

# But I think we need to test a few degrees of heterogeneity in the variation
sd_study_sd<-sd_study_mu * c(0.1, 0.5, 1)

# We will assume a small, medium and large meta-analysis
k_study<-c(30, 100)

# We will drop the SD from, numbers_lost proportion of the data
proportion_lost<-0

# Number of reps for this array of the simulation - it is run on an array of 1000 cores, therefore 10 reps per core
n_reps<-100

# Get all the parameter combinations to run
parameters<-expand.grid(lnRR=lnRR, k_effect_mu=k_effect_mu, k_effect_sd=k_effect_sd, tau2=tau2, icc_study=icc_study, mu_control=mu_control, n_study_mu=n_study_mu, n_study_sd=n_study_sd, sd_study_mu=sd_study_mu, sd_study_sd=sd_study_sd, k_study=k_study, proportion_lost=proportion_lost)

# Remove the parameters with the combinations of low sd n and high mean sd n etc
drop<-c(which(parameters$n_study_mu == 10 & parameters$n_study_sd == sqrt(1.5 * 30)), which(parameters$n_study_mu == 30 & parameters$n_study_sd == sqrt(1.5 * 10)))
parameters<-parameters[-drop,]
rownames(parameters)<-seq(nrow(parameters))

# Create a template to hold the simulation results. For each of the replications, we will perform the analysis by the three methods and record the est, the SE and the tau, then compare to above
template<-data.frame(conv_ests=rep(NA, n_reps), conv_SE=NA, conv_Tau2=NA, wght_ests=NA, wght_SE=NA, wght_Tau2=NA, lnCV1_Tau2=NA, lnCV1_Q=NA, lnCV1_Q_p=NA, lnCV1_df=NA, lnCV2_Tau2=NA, lnCV2_Q=NA, lnCV2_Q_p=NA, lnCV2_df=NA, MA_ests=NA, MA_SE=NA, MA_Tau2=NA)

# Output list to fill in there will be an object for each set of parameters in the batch - that object will be the filled in template
results<-list()

# Loop for the each parameter set
for(p in 1:nrow(parameters)){

	# Get the pth set of parameters
	parameters_p<-parameters[p,]
				
	# Create a copy of the template to fill in for this set of parameters
	template_p<-template
	
	# We will do it reps times, each time deleting n numbers_lost of the SDs and reanalysing by the 3 methods
	for(n in 1:n_reps){
		
		# Keep doing until we get the whole nth row of template_p filled in - repeats simulation in the event that any model has convergence issues
		while(sum(is.na(template_p[n,])) > 0){
		
			# Create an nth data set
			data_n<-sim_data(lnRR=parameters_p$lnRR, k_study=parameters_p$k_study, k_effect_mu=parameters_p$k_effect_mu, k_effect_sd=parameters_p$k_effect_sd, tau2=parameters_p$tau2, icc_study=parameters_p$icc_study, n_study_mu=parameters_p$n_study_mu, n_study_sd=parameters_p$n_study_sd, sd_study_mu=parameters_p$sd_study_mu, sd_study_sd=parameters_p$sd_study_sd)
			
			# Calculate the CVs and effect sizes for the complete data and analyse
			data_n$Control.CV<-data_n$Control.SD / data_n$Control.Mean
			data_n$Treatment.CV<-data_n$Treatment.SD / data_n$Treatment.Mean
			data_n<-cbind(data_n, my_lnRR(Control.Mean=data_n$Control.Mean, Treatment.Mean=data_n$Treatment.Mean, Control.CV=data_n$Control.CV, Treatment.CV=data_n$Treatment.CV, Control.n=data_n$Control.n, Treatment.n=data_n$Treatment.n))
			
			# If we are dropping data find the row IDs you want to have missing SDs and set SD, CV and effect sizes as NA
			# Also create a version of the dataset which has deleted missing SD for complete cases analysis
			if(parameters_p$proportion_lost > 0){
				drop<-sample(seq(1, nrow(data_n), 1), round(nrow(data_n) * parameters_p$proportion_lost))
				data_n$Control.SD[drop]<-NA
				data_n$Treatment.SD[drop]<-NA
				data_n$Control.CV[drop]<-NA
				data_n$Treatment.CV[drop]<-NA
				data_n$yi[drop]<-NA
				data_n$vi[drop]<-NA
				data_drop<-data_n[-drop,]
			}else{
				data_drop<-data_n
			}
			data_drop<-droplevels(data_drop)
			
			# Start with the conventional analysis - when we have missing data this equates to a complete cases analysis
			conv_model<-try(rma.mv(yi = yi, V = vi, random=list(~1|Effect), data=data_drop), silent=T)
			
			# Save the results if the model fitted
			if(class(conv_model)[1] == "rma.mv"){
				template_p$conv_ests[n]<-conv_model$b
				template_p$conv_SE[n]<-conv_model$se
				template_p$conv_Tau2[n]<-conv_model$sigma2[1]
			}
			
			# Method from the Ecology Letters Paper
			# Use the weighted mean cv everywhere, then proceed to fit a meta-analysis as normal.

			# Here we calculate the weighted mean CV from those data available - this is pooled by study first - note here study = effect - this is a big deal in the MLMA sim, but here is the same as averaging data_n
			pooled<-ddply(data_n, .(Study), summarise, mean_Control.CV=mean(Control.CV, na.rm=T), mean_Treatment.CV=mean(Treatment.CV, na.rm=T), mean_Control.n=mean(Control.n), mean_Treatment.n=mean(Treatment.n)) 			
			mcv1<-weighted.mean(pooled$mean_Treatment.CV, pooled$mean_Treatment.n, na.rm=TRUE)
			mcv2<-weighted.mean(pooled$mean_Control.CV, pooled$mean_Control.n, na.rm=TRUE)
														
			# Now estimate a new vi using the average CVs everywhere - note for consistency yi is used throughout
			data_n$vi_2<-mcv1^2/data_n$Treatment.n + mcv2^2/data_n$Control.n + mcv1^4/(2 * data_n$Treatment.n^2) + mcv2^4/(2 * data_n$Control.n^2)
			if(parameters_p$proportion_lost > 0){
				data_n$yi[drop]<-log(data_n$Treatment.Mean[drop] / data_n$Control.Mean[drop]) + 0.5 * (mcv1/data_n$Treatment.n[drop] + mcv2/data_n$Control.n[drop])
			}
			
			# Drop mcvs to be safe
			rm(mcv1)
			rm(mcv2)
			
			# Fit the model
			wght_model<-try(rma.mv(yi = yi, V = vi_2, random=list(~1|Effect), data=data_n), silent=T)
			
			# Save the results if the model fitted
			if(class(wght_model)[1] == "rma.mv"){
				template_p$wght_ests[n]<-wght_model$b
				template_p$wght_SE[n]<-wght_model$se
				template_p$wght_Tau2[n]<-wght_model$sigma2[1]
			}

			# Now the two-step meta-analysis, starting with analysis of lnCVR - need to calculate the new effect sizes first
			data_lnCV1<-cbind(data_drop[,-c(11:12)], my_lnCV(Mean=data_drop$Treatment.Mean, SD=data_drop$Treatment.SD, n=data_drop$Treatment.n))
			data_lnCV2<-cbind(data_drop[,-c(11:12)], my_lnCV(Mean=data_drop$Control.Mean, SD=data_drop$Control.SD, n=data_drop$Control.n))			

			# Fit the models for lnCV in each group
			lnCV1_model<-try(rma.mv(yi = yi, V = vi, random=list(~1|Effect), data=data_lnCV1), silent=T)
			lnCV2_model<-try(rma.mv(yi = yi, V = vi, random=list(~1|Effect), data=data_lnCV2), silent=T)			

			# If they both fitted save the results and progress to step 2 of the analyses
			if(class(lnCV1_model)[1] == "rma.mv" & class(lnCV2_model)[1] == "rma.mv"){
				
				# Save the results of the lnCV1 model - could be useful
				template_p$lnCV1_Tau2[n]<-lnCV1_model$sigma2[1]
				template_p$lnCV1_Q[n]<-lnCV1_model$QE
				template_p$lnCV1_Q_p[n]<-lnCV1_model$QEp
				template_p$lnCV1_df[n]<-nrow(data_lnCV1)-1

				# Save the results of the lnCV2 model - could be useful
				template_p$lnCV2_Tau2[n]<-lnCV2_model$sigma2[1]
				template_p$lnCV2_Q[n]<-lnCV2_model$QE
				template_p$lnCV2_Q_p[n]<-lnCV2_model$QEp
				template_p$lnCV2_df[n]<-nrow(data_lnCV2)-1
				
				# Calculate new vi's for lnRR using the estimates from lnCV1 and lnCV2 this is where the real action is!
				mcv1<-exp(lnCV1_model$b + 0.5 * lnCV1_model$sigma2[1])[1]
				mcv2<-exp(lnCV2_model$b + 0.5 * lnCV2_model$sigma2[1])[1]
				data_n$vi_3<-mcv1^2/data_n$Treatment.n + mcv2^2/data_n$Control.n + mcv1^4/(2 * data_n$Treatment.n^2) + mcv2^4/(2 * data_n$Control.n^2)
				if(parameters_p$proportion_lost > 0){
					data_n$yi[drop]<-log(data_n$Treatment.Mean[drop] / data_n$Control.Mean[drop]) + 0.5 * (mcv1/data_n$Treatment.n[drop] + mcv2/data_n$Control.n[drop])
				}
				
				# Drop mcvs to be safe
				rm(mcv1)
				rm(mcv2)
	
				# Fit the model for lnRR
				MA_model<-try(rma.mv(yi = yi, V = vi_3, random=list(~1|Effect), data=data_n), silent=T)

				# Save the results if the model fitted
				if(class(MA_model)[1] == "rma.mv"){
					template_p$MA_ests[n]<-MA_model$b
					template_p$MA_SE[n]<-MA_model$se
					template_p$MA_Tau2[n]<-MA_model$sigma2[1]
				}
	
			}# Closing the outer if statement for the lnCV models

		}# Closing the while loop
			
	}# Closing the reps
	
	# Save the results from this n_reps set of parameters
	results[[p]]<-template_p

}# Closing the parameters loop

# Save the results from the simulation
save(results, file=paste0("Simulation_1_", index, ".Rdata"))
save(parameters, file="Parameters_Simulation_1.Rdata")



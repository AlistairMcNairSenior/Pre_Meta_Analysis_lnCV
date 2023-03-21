
# This script was used to aggregate multiple instances from the array job on the HPC. It also then calculates the bias, coverage and aggregate statsitics used in the plotting

# Clean up the R Environment 
rm(list=ls())

# Set the working directory
wd<-"/Users/alistairsenior/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Step1_lnCV/HPC Simulation"
setwd(wd)

# Load the relevant libraries, and header
source("0.Header.R")
library(plyr)

# Check the files
files<-dir("Output_Simulation_1")[-c(1:4)]

load(paste0("Output_Simulation_1/", files[1]))
res_unlist<-results
for(i in 2:length(files)){
	load(paste0("Output_Simulation_1/", files[i]))
	for(j in 1:length(res_unlist)){
		res_unlist[[j]]<-rbind(res_unlist[[j]], results[[j]])
	}
	print(i)
}
results<-res_unlist

# Load the parameters 
load("Output_Simulation_1/Parameters_Simulation_1.Rdata")

# Long format all the results
for(i in 1:length(results)){
	long_i<-cbind(results[[i]], parameters[i,])
	if(i == 1){
		long<-long_i
	}else{
		long<-rbind(long, long_i)	
	}
}
long_mat<-as.matrix(long)


# Pull out the lnCV data
lnCV_long<-as.data.frame(long_mat[,c(7:14,18:29)])
long_mat<-long_mat[,-c(7:14)]
lnRR_long<-long_mat

# Save that
save(lnRR_long, file="long_lnRR_simulation_1.Rdata")
save(lnCV_long, file="long_lnCV_simulation_1.Rdata")

# Reformat from long to wide for the methods
long_full<-as.data.frame(rbind(long_mat[,-c(4:9)], long_mat[,-c(1:3,7:9)], long_mat[,-c(1:6)]))
long_full$Method<-c(rep("Conventional", dim(long_mat)[1]), rep("n Weighted", dim(long_mat)[1]), rep("MA Weighted", dim(long_mat)[1]))
names(long_full)[c(1:3)]<-c("Ests", "SE", "Tau2")

# Calculate the bias and coverage
long_full$bias<-long_full$Est - long_full$lnRR
dfs<-long_full$k_study-1
ts<-qt(0.975, df=dfs)
long_full$coverage_t<-(abs(long_full$bias) - ts*long_full$SE) <= 0
long_full$coverage_z<-(abs(long_full$bias) - 1.96*long_full$SE) <= 0
long_full$bias_tau2<-long_full$Tau2 - long_full$tau2

# Add a unique code for each parameter set
long_full$code<-paste0(long_full[,4], "_", long_full[,5])
lnCV_long$code<-paste0(lnCV_long[,9], "_", lnCV_long[,10])
for(i in 6:15){
	long_full$code<-paste0(long_full$code, "_", long_full[,i])
	lnCV_long$code<-paste0(lnCV_long$code, "_", lnCV_long[,i+5])
	print(i)
}

# Now get the aggregate stats
long_full$code2<-paste0(long_full$code, long_full$Method)
agg_stats<-ddply(long_full, .(code2), summarise, code=code[1], Method=Method[1], lower_bias=my_median(bias)[1], median_bias=my_median(bias)[2], upper_bias=my_median(bias)[3], lower_bias_tau2=my_median(bias_tau2)[1], median_bias_tau2=my_median(bias_tau2)[2], upper_bias_tau2=my_median(bias_tau2)[3], coverage_t=mean(coverage_t), coverage_z=mean(coverage_z)) 

# Combine with the parameter set by code
parameters$code<-paste0(parameters[,1], "_", parameters[,2])
for(i in 3:12){
	parameters$code<-paste0(parameters$code, "_", parameters[,i])
	print(i)
}
agg_stats<-cbind(agg_stats, parameters[match(agg_stats$code, parameters$code),])
agg_results<-agg_stats[,-c(1:2)]

# Save that to the githib repo
save(agg_results, file="agg_results_simulation_1.Rdata")


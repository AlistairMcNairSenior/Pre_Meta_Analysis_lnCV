
# Clean up the R Environment 
rm(list=ls())

# Set the working directory
wd<-"/Users/alistairsenior/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Step1_lnCV/HPC Simulation" # MacMini Home
setwd(wd)

# Load the relevant libraries, and header
library(ggplot2)
library(plyr)
library(gridExtra)
library(ggbeeswarm)
library(ggcorrplot)
library(ggthemes)

# Load the results
load("agg_results_simulation_2.Rdata")

head(agg_results)

# The tau2s and ICCs in the simulation - put each combo on a seperate plot - ns and ks on the plot
tau2s<-c(9e-04, 9e-02)
ns<-c(10, 30)
ks<-c(30, 100)
 
# upper Y Limits
ylnR<-c(0.001, 0.005)
yCov<-c(0.98, 0.98)
yTau<-c(0.003, 0.005)
y_counter<-0
		
# Loop for tau2s
for(h in 1:length(tau2s)){
	
	# Create the file for plotting
	y_counter<-y_counter+1
	file_name<-paste0("simulation_2_", tau2s[h], ".pdf")
	pdf(file_name, height=10, width=10)
	
	par(mfrow=c(4,3), mar=c(2,2,2,2), oma=c(4,5,4,1))
	location<-"topright"
	if(h == 2){
		location<-"bottomright"
	}
	
	for(j in 1:length(ks)){
	
		for(i in 1:length(ns)){
				
			plot<-data.frame(tau2=tau2s[h], n_study_mu=ns[i], k_study=ks[j])
			tag<-which(agg_results$tau2 == plot$tau2 & agg_results$n_study_mu == plot$n_study_mu & agg_results$k_study == plot$k_study)
			plot_data<-agg_results[tag,]
			plot_data<-plot_data[order(plot_data$sd_study_sd, plot_data$Method),]
			plot_data$x<-sort(rep(c(1:3), 3)) + c(-0.2, 0, 0.2)	
			print(plot_data)
			
			plot(-100, -100, xlim=c(0.6, 3.4), ylim=c(-ylnR[y_counter], ylnR[y_counter]), cex=2, xaxt="n", xlab="", ylab="", yaxt="n")
			if(i==2 & j == 2){
				mtext(c("Low", "Medium", "High"), at=c(1:3), side=1, line=0.5)
			}
			if(i==1 & j == 1){
				mtext("Bias lnRR", side=3, line=1, cex=1.8)
				legend(location, c("Conv.", "pre-MA", "n-wght"), pch=16, col=c(1:3), cex=1.3)
				
			}
			abline(h=0, lty=2, lwd=2, col="light grey")
			points(plot_data$x, plot_data$median_bias, pch="-", col=as.factor(plot_data$Method), cex=2)
			arrows(plot_data$x, plot_data$lower_bias, plot_data$x, plot_data$upper_bias, col=as.factor(plot_data$Method), code=0, lwd=2)
			axis(side=2, at=seq(-ylnR[y_counter], ylnR[y_counter], length=3), cex.axis=1.5, labels=as.character(seq(-ylnR[y_counter], ylnR[y_counter], length=3)))
			mtext(paste0("n = ", ns[i]), side=2, cex=1.8, line=5)
			mtext(paste0("K = ", ks[j]), side=2, cex=1.8, line=2.75)
			
			plot(-100, -100, xlim=c(0.6, 3.4), ylim=c(0.95 - (yCov[y_counter]-0.95), yCov[y_counter]), cex=2, xaxt="n", xlab="", ylab="", yaxt="n")
			if(i==2 & j == 2){
				mtext(c("Low", "Medium", "High"), at=c(1:3), side=1, line=0.5)
				mtext("Heterogeneity in SD", side=1, line=3.5, cex=1.8)
			}
			if(i==1 & j == 1){
				mtext("Coverage", side=3, line=1, cex=1.8)
			}
			abline(h=0.95, lty=2, lwd=2, col="light grey")
			points(plot_data$x, plot_data$coverage_z, col=as.factor(plot_data$Method), pch=16, cex=2)
			axis(side=2, at=seq(0.95 - (yCov[y_counter]-0.95), yCov[y_counter], length=3), cex.axis=1.5, labels=as.character(seq(0.95 - (yCov[y_counter]-0.95), yCov[y_counter], length=3)))
			
			plot(-100, -100, xlim=c(0.6, 3.4), ylim=c(-yTau[y_counter], yTau[y_counter]), xaxt="n", xlab="", ylab="", yaxt="n")
			arrows(plot_data$x, plot_data$lower_bias_tau2, plot_data$x, plot_data$upper_bias_tau2, col=as.factor(plot_data$Method), code=0, lwd=2)
			if(i==2 & j == 2){
				mtext(c("Low", "Medium", "High"), at=c(1:3), side=1, line=0.5)
			}
			if(i==1 & j == 1){
				mtext("Bias Tau2", side=3, line=1, cex=1.8)
			}
			abline(h=0, lty=2, lwd=2, col="light grey")
			points(plot_data$x, plot_data$median_bias_tau2, col=as.factor(plot_data$Method), pch="-",  cex=2)
			axis(side=2, at=seq(-yTau[y_counter], yTau[y_counter], length=3), cex.axis=1.5, labels=as.character(seq(-yTau[y_counter], yTau[y_counter], length=3)))
		
		}
	
	}
	
	dev.off()

}


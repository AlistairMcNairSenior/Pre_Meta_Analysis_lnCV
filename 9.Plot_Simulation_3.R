

# This script creates a figure for the complete set of all simlulation parameters explored in the third simulation.

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
load("agg_results_simulation_3.Rdata")
head(agg_results)

# The tau2s and ICCs in the simulation - put each combination ns and k on a seperate panel, with different levels of heterogeneity (tau2) on a seperate pdf file
tau2s<-c(9e-02)
ns<-c(10, 30)
ks<-c(30, 100)
 
# Set the upper limits for the y-axis on each plot - found by trial and error 
ylnR<-c(0.005)
yCov<-c(0.97)
yTau<-c(0.005)
yICC<-c(0.02)
y_counter<-0
		
# Loop that makes a plot for each level of tau2
for(h in 1:length(tau2s)){
	
	# Create the file for plotting
	y_counter<-y_counter+1
	file_name<-paste0("simulation_3_", tau2s[h], ".pdf")
	pdf(file_name, height=10, width=12.5)
	
	# visuals for the plot, and how to layout
	par(mfrow=c(4,4), mar=c(2,2,2,2), oma=c(4,5,4,1))
	location<-"bottomright"
	
	# Loop to make plot for each level of k
	for(j in 1:length(ks)){
	
		# Loop to make plot for each level of n
		for(i in 1:length(ns)){
			
			# Subset the the data for hth tau, jth k and ith n for plotting	
			plot<-data.frame(tau2=tau2s[h], n_study_mu=ns[i], k_study=ks[j])
			tag<-which(agg_results$tau2 == plot$tau2 & agg_results$n_study_mu == plot$n_study_mu & agg_results$k_study == plot$k_study)
			plot_data<-agg_results[tag,]
			plot_data<-plot_data[order(plot_data$sd_study_sd, plot_data$Method),]
			
			# Create an x variable for plotting on eahc panel
			plot_data$x<-sort(rep(c(1:3), 3)) + c(-0.2, 0, 0.2)	
			print(plot_data)
			
			# Make the plot for bias in the overall effect size
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
			
			# Make the plot for coverage
			plot(-100, -100, xlim=c(0.6, 3.4), ylim=c(0.95 - (yCov[y_counter]-0.95), yCov[y_counter]), cex=2, xaxt="n", xlab="", ylab="", yaxt="n")
			if(i==2 & j == 2){
				mtext(c("Low", "Medium", "High"), at=c(1:3), side=1, line=0.5)
				mtext("Heterogeneity in SD", side=1, line=3.5, cex=1.8, at=4)
			}
			if(i==1 & j == 1){
				mtext("Coverage", side=3, line=1, cex=1.8)
			}
			abline(h=0.95, lty=2, lwd=2, col="light grey")
			points(plot_data$x, plot_data$coverage_z, col=as.factor(plot_data$Method), pch=16, cex=2)
			axis(side=2, at=seq(0.95 - (yCov[y_counter]-0.95), yCov[y_counter], length=3), cex.axis=1.5, labels=as.character(seq(0.95 - (yCov[y_counter]-0.95), yCov[y_counter], length=3)))
			
			# Make the plot for bias in the heterogeneity
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
			
			# Make the plot for bias in the ICC
			plot(-100, -100, xlim=c(0.6, 3.4), ylim=c(-yICC[y_counter], yICC[y_counter]), xaxt="n", xlab="", ylab="", yaxt="n")
			arrows(plot_data$x, plot_data$lower_bias_icc, plot_data$x, plot_data$upper_bias_icc, col=as.factor(plot_data$Method), code=0, lwd=2)
			if(i==2 & j == 2){
				mtext(c("Low", "Medium", "High"), at=c(1:3), side=1, line=0.5)
			}
			if(i==1 & j == 1){
				mtext("Bias ICC", side=3, line=1, cex=1.8)
			}
			abline(h=0, lty=2, lwd=2, col="light grey")
			points(plot_data$x, plot_data$median_bias_icc, col=as.factor(plot_data$Method), pch="-",  cex=2)
			axis(side=2, at=seq(-yICC[y_counter], yICC[y_counter], length=3), cex.axis=1.5, labels=as.character(seq(-yICC[y_counter], yICC[y_counter], length=3)))
		
		}
	
	}
	
	# Close the file for plotting the hth tau2
	dev.off()

}

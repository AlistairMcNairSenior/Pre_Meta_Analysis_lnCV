
# Clean Up
rm(list=ls())

# Load packages and functions
library(metafor)
source("0.Header.R")
library(orchaRd)
library(gridExtra)
library(ggbeeswarm)

# Load the body mass data
data<-read.csv("Data_Worked_Example_2.csv")
head(data)

# Calculate the CVs
data$Cont_CV<-data$Cont_SD/data$Cont_Mean
data$Treat_CV<-data$Treat_SD/data$Treat_Mean

# Create a unit-level variable to identify each effect size
data$ES<-paste0("ES_", as.factor(seq(1, nrow(data), 1)))

# Calculate Geary's Index following equation 20 in text		
geary_Control<-1/data$Cont_CV * ((4 * data$Cont_n^(3/2)) / (1 + 4*data$Cont_n))
geary_Treat<-1/data$Treat_CV * ((4 * data$Treat_n^(3/2)) / (1 + 4*data$Treat_n))

# No effect sizes failing Geary's test
which(geary_Control < 3) / nrow(data)
which(geary_Treat < 3) / nrow(data)

# Standard meta-analysis of the lnRR
data_lnRR<-cbind(data, my_lnRR(Control.Mean=data$Cont_Mean, Control.CV=data$Cont_CV, Control.n=data$Cont_n, Treatment.Mean=data$Treat_Mean, Treatment.CV=data$Treat_CV, Treatment.n=data$Treat_n))

# Creating the variance-covariance matrix for the shared controls issue
V <- matrix(0,nrow = dim(data_lnRR)[1],ncol = dim(data_lnRR)[1])
rownames(V) <- data$ES
colnames(V) <- data$ES
Shared_ids<-unique(data_lnRR$Control_id[which(data_lnRR$Shared_control == "yes")])
tag<-match(Shared_ids, data_lnRR$Control_id)
shared_data<-data_lnRR[tag,]
shared_cov<-(shared_data$Cont_CV^2 / shared_data$Cont_n) + (shared_data$Cont_CV^4 / (2*shared_data$Cont_n^2))
	
for(i in 1:length(Shared_ids)){
	tags<-which(data_lnRR$Control_id == Shared_ids[i])
	combinations<-combn(tags, 2)
	for(j in 1:dim(combinations)[2]){
		p1 <- combinations[1,j]
  		p2 <- combinations[2,j]
		V[p1,p2] <- shared_cov[i]
		V[p2,p1] <- shared_cov[i]
	}
}

# Add the sampling variances in to the diag of the var-cov matrix for shared controls
diag(V)<-data_lnRR$vi

# Fit the MLMA using the conventional analysis 
conv<-rma.mv(yi=yi, V=V, random=list(~1|Study, ~1|Strain, ~1|ES), data=data_lnRR)

# Pre-meta-analysis - dropping missing SD at this step
# Start with the pre-analysis of the treatment data
lnCV_T<-cbind(data, my_lnCV(Mean=data$Treat_Mean, SD=data$Treat_SD, n=data$Treat_n))
pre_T<-rma.mv(yi=yi, V=vi, random=list(~1|Study, ~1|Strain, ~1|ES), data=lnCV_T)
pre_T
sum(pre_T$sigma2)

# Calculate the I2, using the typical sampling variance
tsvT<-((nrow(lnCV_T) - 1) * sum(1/lnCV_T$vi)) / (sum(1/lnCV_T$vi)^2 - sum((1/lnCV_T$vi)^2))
sum(pre_T$sigma2) / (sum(pre_T$sigma2) + tsvT) * 100

# Partition I2
pre_T$sigma2 / (sum(pre_T$sigma2) + tsvT) * 100

# Repeat for the control Data
data_cont<-data[match(unique(data$Control_id), data$Control_id),]
lnCV_C<-cbind(data_cont, my_lnCV(Mean=data_cont$Cont_Mean, SD=data_cont $Cont_SD, n=data_cont $Cont_n))
pre_C<-rma.mv(yi=yi, V=vi, random=list(~1|Study, ~1|Strain, ~1|ES), data=lnCV_C)
pre_C
sum(pre_C$sigma2)

# I2 for pre control lnCV
tsvC<-((nrow(lnCV_C) - 1) * sum(1/lnCV_C$vi)) / (sum(1/lnCV_C$vi)^2 - sum((1/lnCV_C$vi)^2))
sum(pre_C$sigma2) / (sum(pre_C$sigma2) + tsvC) * 100

# Partition I2
pre_C$sigma2 / (sum(pre_C$sigma2) + tsvC) * 100

# Get the estimates of the CV in the two groups from the pre-meta-analysis
CV_T<-exp(pre_T$b[1] + 0.5*sum(pre_T$sigma2))
CV_C<-exp(pre_C$b[1] + 0.5*sum(pre_C$sigma2))

# Calculate the sampling variances using the new estimates of CV from the pre-meta-analysis
data_lnRR$vi_2<-CV_T^2/data$Treat_n + CV_C^2/data$Cont_n + CV_T^4/(2 * data$Treat_n^2) + CV_C^2/(2 * data$Cont_n^2)

# re-create the vcv for the shared controls
V <- matrix(0,nrow = dim(data_lnRR)[1],ncol = dim(data_lnRR)[1])
rownames(V) <- data$ES
colnames(V) <- data$ES

Shared_ids<-unique(data_lnRR$Control_id[which(data_lnRR$Shared_control == "yes")])
tag<-match(Shared_ids, data_lnRR$Control_id)
shared_data<-data_lnRR[tag,]
shared_cov<-(CV_C^2 / shared_data$Cont_n) + (CV_C^4 / (2*shared_data$Cont_n^2))
	
for(i in 1:length(Shared_ids)){
	tags<-which(data_lnRR$Control_id == Shared_ids[i])
	combinations<-combn(tags, 2)
	for(j in 1:dim(combinations)[2]){
		p1 <- combinations[1,j]
  		p2 <- combinations[2,j]
		V[p1,p2] <- shared_cov[i]
		V[p2,p1] <- shared_cov[i]
	}
}

# Add the new sampling variance in to the diag for var-cov matrix
diag(V)<-data_lnRR$vi_2

# Refit with the MLMA for the lnRR with the new sampling variances
pre_MA<-rma.mv(yi=yi, V=V, random=list(~1|Study, ~1|Strain, ~1|ES), data=data_lnRR)

# Compare the results
conv
pre_MA

# Total heterogeneity
sum(conv$sigma2)
sum(pre_MA$sigma2)

# Calculate I2s
tsv1<-((nrow(data_lnRR) - 1) * sum(1/data_lnRR$vi)) / (sum(1/data_lnRR$vi)^2 - sum((1/data_lnRR$vi)^2))
tsv2<-((nrow(data_lnRR) - 1) * sum(1/data_lnRR$vi_2)) / (sum(1/data_lnRR$vi_2)^2 - sum((1/data_lnRR$vi_2)^2))

# Total I2
round(sum(conv$sigma2) / (sum(conv$sigma2) + tsv1) * 100)
round(sum(pre_MA$sigma2) / (sum(pre_MA$sigma2) + tsv2) * 100)

# Partitioned I2
round(conv$sigma2 / (sum(conv$sigma2) + tsv1) * 100)
round(pre_MA$sigma2 / (sum(pre_MA$sigma2) + tsv2) * 100)

# Make figure 3 orchanrd plots
data_lnRR$N<-data_lnRR$Cont_n + data_lnRR$Treat_n
data_lnRR$lnRR<-data_lnRR$yi
results_conv<-mod_results(conv, mod="1", at=NULL, data=data_lnRR, group="Study")
results_conv
results_pre<-mod_results(pre_MA, mod="1", at=NULL, data=data_lnRR, group="Study")
results_pre

p1<-ggplot(data = data_lnRR, aes(x = lnRR, y = 1, size = N)) +
	geom_vline(xintercept=0) + 
	geom_quasirandom(groupOnX=F, alpha=0.5) + theme_bw() +	
	theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
	axis.title.x=element_text(size=15), axis.text.x=element_text(size=15)) + 
	geom_point(aes(x=results_conv$mod_table$estimate, y=1.05), inherit.aes=F, size=11, colour="gold2", alpha=0.05) +
	geom_point(aes(x=results_pre $mod_table$estimate, y=0.95), inherit.aes=F, size=11, colour="cornflowerblue", alpha=0.05) +
	geom_linerange(aes(xmin=results_conv$mod_table$lowerCL, xmax=results_conv$mod_table$upperCL, y=1.05), inherit.aes=F, colour="gold2", linewidth=3) +
	geom_linerange(aes(xmin=results_conv$mod_table$lowerPR, xmax=results_conv$mod_table$upperPR, y=1.05), inherit.aes=F, colour="gold2", linewidth=0.5) +
	geom_linerange(aes(xmin=results_pre$mod_table$lowerCL, xmax=results_pre$mod_table$upperCL, y=0.95), inherit.aes=F, colour="cornflowerblue", linewidth=3) +
	geom_linerange(aes(xmin=results_pre$mod_table$lowerPR, xmax=results_pre$mod_table$upperPR, y=0.95), inherit.aes=F, colour="cornflowerblue", linewidth=0.5) +
	theme(legend.position=c(0.05, 0.82), legend.background=element_blank(), legend.text=element_text(size=15), legend.title=element_text(size=15))


pdf("Orchard_Plot.pdf", height=5, width=11)
grid.arrange(p1, nrow=1)
dev.off()

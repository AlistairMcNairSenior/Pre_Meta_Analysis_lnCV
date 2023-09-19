
# Clean Up
rm(list=ls())

# Load packages and functions
library(metafor)
source("0.Header.R")

# Load the biodiversity data
data<-read.csv("Data_Worked_Example_3.csv")
head(data)

# Study 20 has 0 SD lets dump
data<-data[-20,]
	
# There are also missing data so lets dump as that is not what we are trying to illustrate here
data<-data[-which(is.na(data$Cont_SD) == T),]
	
# Calculate the CVs
data$Cont_CV<-data$Cont_SD/data$Cont_Mean
data$Treat_CV<-data$Treat_SD/data$Treat_Mean
data$ES<-as.factor(seq(1, nrow(data), 1))

# Calculate Geary's Index following equation 20 in text			
geary_Control<-1/data$Cont_CV * ((4 * data$Cont_n^(3/2)) / (1 + 4*data$Cont_n))
geary_Treat<-1/data$Treat_CV * ((4 * data$Treat_n^(3/2)) / (1 + 4*data$Treat_n))

# Check those failing
length(which(geary_Control < 3)) / nrow(data)
length(which(geary_Treat < 3)) / nrow(data)

# Find those failing
tag<-unique(c(which(geary_Control < 3), which(geary_Treat < 3)))
length(tag) / nrow(data)

# Standard meta-analysis of the lnRR
data_full<-cbind(data, my_lnRR(Control.Mean=data$Cont_Mean, Control.CV=data$Cont_CV, Control.n=data$Cont_n, Treatment.Mean=data$Treat_Mean, Treatment.CV=data$Treat_CV, Treatment.n=data$Treat_n))
# recalcualte yi as the raw lnRR
data_full$yi<-log(data_full$Treat_Mean / data_full$Cont_Mean)

# Subset for those passing
data_pass<-data_full[-tag,]

# Fit model to the full data and to those which pass Geary
conv<-rma.mv(yi=yi, V=vi, random=list(~1|ES, ~1|study.name), data=data_full)
conv_sensi<-rma.mv(yi=yi, V=vi, random=list(~1|ES, ~1|study.name), data=data_pass)

# Pre-meta-analysis - dropping missing SD at this step
# Start with the pre-analysis of the treatment data
lnCV_T<-cbind(data_full[ ,-c(36,37)], my_lnCV(Mean=data_full$Treat_Mean, SD=data_full$Treat_SD, n=data_full$Treat_n))
pre_T<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_T)

# Now do the pre-meta-analysis for the control data
lnCV_C<-cbind(data_full[ ,-c(36,37)], my_lnCV(Mean=data_full$Cont_Mean, SD=data_full$Cont_SD, n=data_full$Cont_n))
pre_C<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_C)

# Check the results
pre_T
pre_C

# Get the estimates of the CV in the two groups from the pre-meta-analysis
CV_T<-exp(pre_T$b[1] + 0.5*sum(pre_T$sigma2))
CV_C<-exp(pre_C$b[1] + 0.5*sum(pre_C$sigma2))

# Calculate the sampling variances using the new estimates of CV from the pre-meta-analysis
data_full$vi_2<-CV_T^2/data_full$Treat_n + CV_C^2/data_full$Cont_n + CV_T^4/(2 * data_full$Treat_n^2) + CV_C^2/(2 * data_full$Cont_n^2)

# Refit with the REMA for the lnRR with the new sampling variances
pre_MA<-rma.mv(yi=yi, V=vi_2, random=list(~1|ES, ~1|study.name), data=data_full)

# Redo the pre-meta-analysis for the data passing Geary
lnCV_T<-cbind(data_pass[ ,-c(36,37)], my_lnCV(Mean=data_pass$Treat_Mean, SD=data_pass$Treat_SD, n=data_pass$Treat_n))
pre_T<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_T)
lnCV_C<-cbind(data_pass[ ,-c(36,37)], my_lnCV(Mean=data_pass$Cont_Mean, SD=data_pass$Cont_SD, n=data_pass$Cont_n))
pre_C<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_C)

# Check the results
pre_T
pre_C

# Get the estimates of the CV in the two groups from the pre-meta-analysis
CV_T<-exp(pre_T$b[1] + 0.5*sum(pre_T$sigma2))
CV_C<-exp(pre_C$b[1] + 0.5*sum(pre_C$sigma2))

# Calculate the variances
data_pass$vi_2<-CV_T^2/data_pass$Treat_n + CV_C^2/data_pass$Cont_n + CV_T^4/(2 * data_pass$Treat_n^2) + CV_C^2/(2 * data_pass$Cont_n^2)

# Calculate the sampling variances using the new estimates of CV from the pre-meta-analysis with those data passing Geary
pre_MA_sensi<-rma.mv(yi=yi, V=vi_2, random=list(~1|ES, ~1|study.name), data=data_pass)

# Compare the results
conv
conv_sensi
pre_MA
pre_MA_sensi



# Clean Up
rm(list=ls())

# Load packages and functions
library(metafor)
source("0.Header.R")

# Load the FBG data
data<-read.csv("Data_Worked_Example_3.csv")
head(data)

# Study 20 has 0 SD lets dump
data<-data[-20,]
	
# There are also missing data so lets dump as that is nto what we are trying to illustrate here
data<-data[-which(is.na(data$Cont_SD) == T),]
	
# Calculate the CVs
data$Cont_CV<-data$Cont_SD/data$Cont_Mean
data$Treat_CV<-data$Treat_SD/data$Treat_Mean
data$ES<-as.factor(seq(1, nrow(data), 1))

# Calculate Geary's Index		
geary_Control<-1/data$Cont_CV * ((4 * data$Cont_n^(3/2)) / (1 + 4*data$Cont_n))
geary_Treat<-1/data$Treat_CV * ((4 * data$Treat_n^(3/2)) / (1 + 4*data$Treat_n))

# No effect sizes failing Geary's test
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

# Fit model to the full data and to that which passes Geary
conv<-rma.mv(yi=yi, V=vi, random=list(~1|ES, ~1|study.name), data=data_full)
conv_sensi<-rma.mv(yi=yi, V=vi, random=list(~1|ES, ~1|study.name), data=data_pass)

# Pre-meta-analysis for the full data
lnCV_T<-cbind(data_full[ ,-c(36,37)], my_lnCV(Mean=data_full$Treat_Mean, SD=data_full$Treat_SD, n=data_full$Treat_n))
pre_T<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_T)

lnCV_C<-cbind(data_full[ ,-c(36,37)], my_lnCV(Mean=data_full$Cont_Mean, SD=data_full$Cont_SD, n=data_full$Cont_n))
pre_C<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_C)

# Check the results
pre_T
pre_C

# Get the estimates
CV_T<-exp(pre_T$b[1] + 0.5*sum(pre_T$sigma2))
CV_C<-exp(pre_C$b[1] + 0.5*sum(pre_C$sigma2))

# Calculate the variances
data_full$vi_2<-CV_T^2/data_full$Treat_n + CV_C^2/data_full$Cont_n + CV_T^4/(2 * data_full$Treat_n^2) + CV_C^2/(2 * data_full$Cont_n^2)

# Refit with the new sampling variances
pre_MA<-rma.mv(yi=yi, V=vi_2, random=list(~1|ES, ~1|study.name), data=data_full)

# Redo the pre-meta-analysis for the passing data
lnCV_T<-cbind(data_pass[ ,-c(36,37)], my_lnCV(Mean=data_pass$Treat_Mean, SD=data_pass$Treat_SD, n=data_pass$Treat_n))
pre_T<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_T)

lnCV_C<-cbind(data_pass[ ,-c(36,37)], my_lnCV(Mean=data_pass$Cont_Mean, SD=data_pass$Cont_SD, n=data_pass$Cont_n))
pre_C<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_C)

# Check the results
pre_T
pre_C

# Get the estimates
CV_T<-exp(pre_T$b[1] + 0.5*sum(pre_T$sigma2))
CV_C<-exp(pre_C$b[1] + 0.5*sum(pre_C$sigma2))

# Calculate the variances
data_pass$vi_2<-CV_T^2/data_pass$Treat_n + CV_C^2/data_pass$Cont_n + CV_T^4/(2 * data_pass$Treat_n^2) + CV_C^2/(2 * data_pass$Cont_n^2)

# Refit with the new sampling variances
pre_MA_sensi<-rma.mv(yi=yi, V=vi_2, random=list(~1|ES, ~1|study.name), data=data_pass)

# Compare the results
conv
conv_sensi
pre_MA
pre_MA_sensi

# Try out the conversion to fix some effect sizes
tag<-which((geary_Control < 3 | geary_Treat < 3) & data$Treat_Mean > 1 & data$Cont_Mean > 1)
ln_Treat_Mean<-log(data$Treat_Mean[tag])
ln_Cont_Mean<-log(data$Cont_Mean[tag])
ln_Treat_SD<-data$Treat_SD[tag] / data$Treat_Mean[tag]
ln_Cont_SD<-data$Cont_SD[tag] / data$Cont_Mean[tag]

# Now add those back in and re-calculate Geary
data$Treat_Mean[tag]<-ln_Treat_Mean
data$Cont_Mean[tag]<-ln_Cont_Mean
data$Treat_SD[tag]<-ln_Treat_SD
data$Cont_SD[tag]<-ln_Cont_SD

# Recalculate CVs
data$Cont_CV<-data$Cont_SD/data$Cont_Mean
data$Treat_CV<-data$Treat_SD/data$Treat_Mean

# Calculate Geary's Index		
geary_Control<-1/data$Cont_CV * ((4 * data$Cont_n^(3/2)) / (1 + 4*data$Cont_n))
geary_Treat<-1/data$Treat_CV * ((4 * data$Treat_n^(3/2)) / (1 + 4*data$Treat_n))

# Find those failing
tag<-unique(c(which(geary_Control < 3), which(geary_Treat < 3)))
length(tag) / nrow(data)
# Only saves a few effect sizes 11% to 10%

# Clean Up
rm(list=ls())

# Load packages and functions
library(metafor)
source("0.Header.R")

# Load the Fasting blood glucose data
data<-read.csv("Data_Worked_Example_1.csv")
head(data)
	
# Label the data by author-year for plotting
data$label<-paste0(data$Author, data$Year)

# Calculate the CVs for each group
data$Cont_CV<-data$Cont_SD/data$Cont_Mean
data$Treat_CV<-data$Treat_SD/data$Treat_Mean
data$ES<-as.factor(seq(1, nrow(data), 1))

# Calculate Geary's Index following equation 20 in text	
geary_Control<-1/data$Cont_CV * ((4 * data$Cont_n^(3/2)) / (1 + 4*data$Cont_n))
geary_Treat<-1/data$Treat_CV * ((4 * data$Treat_n^(3/2)) / (1 + 4*data$Treat_n))

# No effect sizes failing Geary's test
which(geary_Control < 3)
which(geary_Treat < 3)

# Standard meta-analysis of the lnRR
data_lnRR<-cbind(data, my_lnRR(Control.Mean=data$Cont_Mean, Control.CV=data$Cont_CV, Control.n=data$Cont_n, Treatment.Mean=data$Treat_Mean, Treatment.CV=data$Treat_CV, Treatment.n=data$Treat_n))

# Sort by effect size
data_lnRR<-data_lnRR[order(data_lnRR$yi),]

# Note there is one with missing SD data
missing_data<-which(is.na(data_lnRR$Cont_SD)==T)

# Fit REMA using the standard method
conv<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=data_lnRR, slab=label)

# Pre-meta-analysis - dropping missing SD at this step
# Start with the pre-analysis of the treatment data
lnCV_T<-cbind(data_lnRR[-missing_data,-c(20,21)], my_lnCV(Mean=data_lnRR$Treat_Mean[-missing_data], SD=data_lnRR$Treat_SD[-missing_data], n=data_lnRR $Treat_n[-missing_data]))
pre_T<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_T)

# Now do the pre-meta-analysis for the control data
lnCV_C<-cbind(data_lnRR[-missing_data,-c(20,21)], my_lnCV(Mean=data_lnRR$Cont_Mean[-missing_data], SD=data_lnRR$Cont_SD[-missing_data], n=data_lnRR$Cont_n[-missing_data]))
pre_C<-rma.mv(yi=yi, V=vi, random=list(~1|ES), data=lnCV_C)

# Check the results
pre_T
pre_C

# I2 for pre treatment lnCV
tsvT<-((nrow(lnCV_T) - 1) * sum(1/lnCV_T$vi)) / (sum(1/lnCV_T$vi)^2 - sum((1/lnCV_T$vi)^2))
pre_T$sigma2 / (pre_T$sigma2 + tsvT) * 100

# I2 for pre control lnCV
tsvC<-((nrow(lnCV_C) - 1) * sum(1/lnCV_C$vi)) / (sum(1/lnCV_C$vi)^2 - sum((1/lnCV_C$vi)^2))
pre_C$sigma2 / (pre_C$sigma2 + tsvC) * 100

# Get the estimates of the CV in the two groups from the pre-meta-analysis
CV_T<-exp(pre_T$b[1] + 0.5*sum(pre_T$sigma2))
CV_C<-exp(pre_C$b[1] + 0.5*sum(pre_C$sigma2))

# Calculate the sampling variances using the new estimates of CV from the pre-meta-analysis
data_lnRR$vi_2<-CV_T^2/data$Treat_n + CV_C^2/data$Cont_n + CV_T^4/(2 * data$Treat_n^2) + CV_C^2/(2 * data$Cont_n^2)
tag<-which(is.na(data_lnRR$yi) == T)
data_lnRR$yi[tag]<-log(data_lnRR$Treat_Mean[tag] / data_lnRR$Cont_Mean[tag])

# Refit with the REMA for the lnRR with the new sampling variances
pre_MA<-rma.mv(yi=yi, V=vi_2, random=list(~1|ES), data=data_lnRR, slab=label)

# Check the results
conv
pre_MA

# Calculate I2s - note inclusion of missing data for conventional analysis
tsv1<-((nrow(data_lnRR[-missing_data]) - 1) * sum(1/data_lnRR$vi[-missing_data])) / (sum(1/data_lnRR$vi[-missing_data])^2 - sum((1/data_lnRR$vi[-missing_data])^2))
tsv2<-((nrow(data_lnRR) - 1) * sum(1/data_lnRR$vi_2)) / (sum(1/data_lnRR$vi_2)^2 - sum((1/data_lnRR$vi_2)^2))

# Total I2
round(sum(conv$sigma2) / (sum(conv$sigma2) + tsv1) * 100)
round(sum(pre_MA$sigma2) / (sum(pre_MA$sigma2) + tsv2) * 100)


# Make some forest plots
pdf("Forest_plots.pdf", height=10, width=13)
par(mfrow=c(1,2))
forest(conv, header="Study", showweights=T, at=seq(-0.4, 0.4, 0.4), xlab="lnRR", cex.lab=1.25)
text(-0.75, -1, "Tau2 = 0.015")
text(0, 21, "Conventional Analysis", font=2, cex=1.5)
forest(pre_MA, header="Study", showweights=T, at=seq(-0.4, 0.4, 0.4), xlab="lnRR", cex.lab=1.25)
text(-0.9, -1, "Tau2 = 0.017")
text(0, 22, "Pre-Meta-Analysis", font=2, cex=1.5)
dev.off()


# In the supplementary I suggest arms-based meta-analysis/meta-regression as another way to do the pre-analysis - this is how to implement that method
# Format the data to longformat for treatment and control groups
lnCV_T$group<-as.factor("1")
lnCV_C$group<-as.factor("0")
lnCV_long<-rbind(lnCV_C, lnCV_T)
lnCV_long$unit<-as.factor(seq(1, nrow(lnCV_long), 1))

# Arm-based meta-analysis, or meta-regression
meta_reg<-rma.mv(yi=yi, V=vi, mods=~group, random=list(~group|ES, ~1|unit), data=lnCV_long)
meta_reg

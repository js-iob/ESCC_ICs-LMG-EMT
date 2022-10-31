#Author: K.T.Shreya Parthasarathi
#Date: 20/09/2022
#Purpose: 5 year survival analysis of patients with ESCC using tcga esophageal cancer ion channel mRNA expression profiles

setwd("path\\to\\current\\working\\directory")

#Importing libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("RTCGA.clinical")

rna <- read.table('ICs_escc_042022.txt', header=T,sep='\t', fill = TRUE)
rna_t = t(rna)
rna_transposed = write.table(rna_t, file="ICs_transposed.txt")

library(RTCGA.clinical)
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)

#Data processing
dim(ESCA.clinical)
rna_ic <- read.table('ICs_transposed.txt', header=T,sep='\t', fill = TRUE)
dim(rna_ic)
rna_ic = distinct(rna_ic)
dim(rna_ic)
colnames(rna_ic)[1] = 'bcr_patient_barcode'
head(rna_ic)

clin = survivalTCGA(ESCA.clinical)
head(clin)
table(clin$patient.vital_status)
#write.table(clin, 'clinical.tsv',sep = '\t')


#Combining gene expression and clinical data
clin = rna_ic %>%
  as_tibble() %>%
  select(bcr_patient_barcode, GJA1, GABRE, GABRR2, GABRQ, GABRA3, ANO1, TRPV3, KCNN4, TRPM7, TRPC1, 
         TRPM2, ITPR3, GABRP, GABRA4) %>%
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
  inner_join(clin, by="bcr_patient_barcode")
names(clin)

dim(clin)
head(clin)
table(clin$patient.vital_status)

#5yr criteria:
all_clin = within(clin, patient.vital_status[patient.vital_status == 1 & times > 1825] <- 0)
all_clin = within(clin, times[times > 1825] <- 1825)
table(all_clin$patient.vital_status)
range (all_clin$times)


#median(all_clin$CLIC1)
#hist(clin$CLIC2)
#median(all_clin$GJA1,na.rm = TRUE)
#Grouping
gene = cut(all_clin$ITPR3, breaks = c(0, median(all_clin$ITPR3,na.rm = TRUE), Inf), labels = c("low", "high"))
sfit = survfit(Surv(times, patient.vital_status)~gene, data = all_clin)
sfit

#Cox regression
fit = coxph(Surv(times, patient.vital_status)~ITPR3, data = all_clin)
fit

#Plotting
ggsurv = ggsurvplot(sfit, legend.labs=c("Low","High"), legend.title=
                      'Expression', title = 'ITPR3', pval = TRUE, pval.method = TRUE, xlab = "Time(in days)" )

ggsurv$plot + ggplot2::annotate("text",x = Inf, y = Inf, vjust = 1, hjust=1, 
                                label = "HR = 4.09 \n p(HR) = 0.04")




#------ Extension anlysis of the project titled "Diversity of Vibrio cholerae O1 through the human gastrointestinal tract during cholera"--------#
#---
#title: "Stool vs Vomit RNA-seq Report of SMIC Study"
#author: - name: "Md Naiem Hossain"
#affiliation: "Research Officer, icddr,b"

#---
  
# Set the working directory

setwd("D:/ICDDRB_Feb/1.RNAseq/1.Final_VC_RNAseq_10_11")
getwd()
#------------------------------------------------------------------
# Load the required packages 
library(limma)
library(edgeR)
library(dplyr)
#------------------------------------------------------------------
# Remooves existing file before proceed 
rm(list=ls());gc()
#------------------------------------------------------------------
# Step 1: Load the count data
count_data_all <- read.csv("ASVS.counts.genes.csv",check.names = F)  # Gene expression matrix
head(count_data_all)
dim(count_data_all)
#-----------------------------------------------------------------
# Step 2: Merge the  count data with annotation file
annotation_file <- read.csv("VC_genes-Ng-file.csv",check.names = F)  # annotation_file
colnames(annotation_file)
head(annotation_file)
dim(annotation_file)
#-----------------------------------------------------------------
#Step 3: Merge the files
Merged_Data <- count_data_all %>%inner_join(annotation_file, by = "Geneid") %>%mutate(GENEID =`N16961 locus`)
dim(Merged_Data)
n_distinct(Merged_Data$GENEID)
n_distinct(Merged_Data$Geneid)
#dim(Merged_Data)
#colnames(Merged_Data)
#-----------------------------------------------------------------
#Step 4:  write the Merged_Data file

write.csv(Merged_Data,"Merged_Data.csv", row.names = F)
#-----------------------------------------------------------------
#Step 5:  read the Merged_Data file ( if required)
#Merged_Data <- read.csv("Merged_Data.csv",check.names = F)
#-----------------------------------------------------------------
# ------------------- Paired Sample Analysis -----------------------
# -----Create output directory For Paired Samples Analysis-------
dir.create("1.Final_Paired_results", showWarnings = FALSE)
# Step P01 : Select the sample names of paired ID along with GENEID
#Paired Sample list 
#Stool -S21,S23,S24,S25,S26
#Vomit -S30,S33,S34,S35,S37
Merged_Data_paired_count <- Merged_Data %>%dplyr::select(GENEID,S21,S23,S24,S25,S26,S30,S33,S34,S35,S37) %>%distinct(GENEID, .keep_all = TRUE)
dim(Merged_Data_paired_count)
write.csv(Merged_Data_paired_count,"1.Final_Paired_results/Merged_Data_paired_count.csv", row.names = F)

#head(Merged_Data_paired_count)
#nrow(Merged_Data_paired_count)
#n_distinct(Merged_Data_paired_count$GENEID)
#nrow(Merged_Data)
#n_distinct(Merged_Data$Geneid)
#Step P02 : Preparing the paired sample list column
Paired_Sample <- data.frame(
  Sample_type = c("Stool", "Stool", "Stool", "Stool", "Stool",
                  "Vomit", "Vomit", "Vomit", "Vomit", "Vomit"),
  Sample = c("S21", "S23", "S24", "S25", "S26",
             "S30", "S33", "S34", "S35", "S37"))

Paired_Sample

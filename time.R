########### time
# keletal muscle
{
# GSE25941_series_matrix.txt
# Step 1: Read the entire file
file_content <- readLines("GSE25941_series_matrix.txt")

# Step 2: Extract metadata (Sample Characteristics)
sample_characteristics_ch2_line <- grep("!Sample_characteristics_ch2", file_content)
sample_characteristics_ch2 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch2_line]), "\t"))[-1]

sample_characteristics_ch3_line <- grep("!Sample_characteristics_ch3", file_content)
sample_characteristics_ch3 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch3_line]), "\t"))[-1]

# Step 3: Read the expression matrix
matrix_start <- grep("!series_matrix_table_begin", file_content)
expression_matrix <- read.delim(textConnection(file_content[(matrix_start + 1):length(file_content)]), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 4: Convert sample characteristics to rows
# Create a new data frame with sample characteristics as rows
metadata_df <- data.frame(
  ID_REF = c("Sample_Characteristics_Ch2", "Sample_Characteristics_Ch3"),
  rbind(sample_characteristics_ch2, sample_characteristics_ch3),
  stringsAsFactors = FALSE
)

# Ensure that the column names of the metadata match the matrix data (excluding ID_REF)
colnames(metadata_df) <- colnames(expression_matrix)

# Step 5: Combine the metadata and the expression matrix
combined_GSE25941 <- rbind(metadata_df, expression_matrix)

# Step 6: Create matrices for each group (young females, young males, old females, old males)
# Exclude the last row using -nrow(combined_GSE25941)
young_females_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Female" & combined_GSE25941[1, ] == "age: Young"])
young_males_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Male" & combined_GSE25941[1, ] == "age: Young"])
old_females_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Female" & combined_GSE25941[1, ] == "age: Old"])
old_males_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Male" & combined_GSE25941[1, ] == "age: Old"])

young_females_matrix <- apply(young_females_matrix, 2, function(x) as.numeric(as.character(x)))
young_males_matrix <- apply(young_males_matrix, 2, function(x) as.numeric(as.character(x)))
old_females_matrix <- apply(old_females_matrix, 2, function(x) as.numeric(as.character(x)))
old_males_matrix <- apply(old_males_matrix, 2, function(x) as.numeric(as.character(x)))

# Step 7: Standardize (scale) the matrices for each group
# This ensures that each feature has a mean of 0 and a standard deviation of 1
young_females_matrix <- scale(young_females_matrix)
young_males_matrix <- scale(young_males_matrix)
old_females_matrix <- scale(old_females_matrix)
old_males_matrix <- scale(old_males_matrix)

# Step 8: Output the dimensions of each matrix to verify correctness
print(dim(young_females_matrix))
print(dim(young_males_matrix))
print(dim(old_females_matrix))
print(dim(old_males_matrix))

# Step 9: Combine young females and young males matrices, and old females and old males matrices
young_combined_matrix <- cbind(young_females_matrix, young_males_matrix)
old_combined_matrix <- cbind(old_females_matrix, old_males_matrix)

# Step 10: Output the dimensions of each combined matrix to verify correctness
print(dim(young_combined_matrix))
print(dim(old_combined_matrix))

group1 <- t(young_combined_matrix)
group2 <- t(old_combined_matrix)

dim(group1)
dim(group2)

library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# SARS-CoV-2
{
# GSE156063_series_matrix.txt
# GSE156063_series_matrix.txt Processing Script
# Step 1: Read the series matrix file (metadata file)
file_content <- readLines("GSE156063_series_matrix.txt")

# Step 2: Extract Sample Characteristics (Sample Characteristics ch2: disease states)
sample_characteristics_ch2_line <- grep("!Sample_characteristics_ch2", file_content)

# Split the characteristics data into individual elements by tab separator
sample_characteristics_ch2 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch2_line]), "\t"))[-1]

# Step 3: Read the gene expression data from the CSV file (GSE156063_swab_gene_counts.csv)
gene_expression_data <- read.csv("GSE156063_swab_gene_counts.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Step 4: Convert sample characteristics to rows
# Create a new data frame with sample characteristics ch2 as a row
metadata_df <- data.frame(
  ID_REF = "Sample_Characteristics_Ch2",  # Label for the metadata row
  t(sample_characteristics_ch2),  # Transpose the sample characteristics into a row format
  stringsAsFactors = FALSE
)

# Step 5: Ensure that the column names of the metadata match the gene expression data
colnames(metadata_df) <- colnames(gene_expression_data)

# Step 6: Combine the metadata and the gene expression data
combined_GSE156063 <- rbind(metadata_df, gene_expression_data)

# Step 7: Save the combined data as a CSV file
# write.csv(combined_GSE156063, "GSE156063_combined_data.csv", row.names = FALSE)
dim(combined_GSE156063)

# Step 8: Create matrices for each disease state group (SC2, other viral, no virus)
# Select samples corresponding to each group based on disease state
sc2_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: SC2"])
other_virus_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: other virus"])
no_virus_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: no virus"])

# Convert character data to numeric if necessary
sc2_matrix <- apply(sc2_matrix, 2, as.numeric)
other_virus_matrix <- apply(other_virus_matrix, 2, as.numeric)
no_virus_matrix <- apply(no_virus_matrix, 2, as.numeric)

# Step 9: Output the dimensions of each matrix to verify
print(dim(sc2_matrix))
print(dim(other_virus_matrix))
print(dim(no_virus_matrix))

# Step 10: Combine SC2 and other viral matrices into a new matrix
sc2_other_combined_matrix <- cbind(sc2_matrix, other_virus_matrix)

# Step 11: Output the dimensions of the combined matrix and the no virus matrix to verify
print(dim(sc2_other_combined_matrix))
print(dim(no_virus_matrix))

group1 <- t(sc2_other_combined_matrix)
group2 <- t(no_virus_matrix)

dim(group1)
dim(group2)

library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# COVID-19 from HDNRA
{
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

data("COVID19")
group1 <- as.matrix(COVID19[c(2:19, 82:87), ])
dim(group1)
group2 <- as.matrix(COVID19[-c(1:19, 82:87), ])
dim(group2)

dim(group1)
dim(group2)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# yeoh-2002-v2
{
# Load the dataset and skip the first row
data1 <- read.table("yeoh-2002-v2_database.txt", header=TRUE, sep="\t", fill=TRUE, skip=1)

# Verify the data structure
# head(data1)
# colnames(data1)

# Create matrices for each sample type
BCR_ABL <- as.matrix(data1[, grep("BCR.ABL", names(data1))])
E2A_PBX1 <- as.matrix(data1[, grep("E2A.PBX1", names(data1))])
Hyperdiploid_50 <- as.matrix(data1[, grep("Hyperdip.50", names(data1))])
MLL <- as.matrix(data1[, grep("MLL", names(data1))])
TEL_AML1 <- as.matrix(data1[, grep("TEL.AML", names(data1))])
T_ALL <- as.matrix(data1[, grep("T.ALL", names(data1))])

# Output the dimensions of each matrix to verify
print(dim(BCR_ABL))
print(dim(E2A_PBX1))
print(dim(Hyperdiploid_50))
print(dim(MLL))
print(dim(TEL_AML1))
print(dim(T_ALL))


group1 <- t(BCR_ABL)
group2 <- t(TEL_AML1)

dim(group1)
dim(group2)

library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# alizadeh-2000-v3
{
# Load the dataset
data2 <- read.table("alizadeh-2000-v3_database.txt", header=TRUE, sep="\t", fill=TRUE, skip=1)

# Verify the data structure
# head(data2)

# Create matrices for each sample type
DLBCL1 <- as.matrix(data2[, grep("DLBCL1", names(data2))])
DLBCL2 <- as.matrix(data2[, grep("DLBCL2", names(data2))])
FL <- as.matrix(data2[, grep("FL", names(data2))])
CLL <- as.matrix(data2[, grep("CLL", names(data2))])

# Convert the matrices to numeric, handling any non-numeric values
DLBCL1 <- apply(DLBCL1, 2, as.numeric)
DLBCL2 <- apply(DLBCL2, 2, as.numeric)
FL <- apply(FL, 2, as.numeric)
CLL <- apply(CLL, 2, as.numeric)

# Convert the matrices to numeric and replace NA values with 0
DLBCL1 <- matrix(ifelse(is.na(as.numeric(DLBCL1)), 0, as.numeric(DLBCL1)), nrow = nrow(DLBCL1), ncol = ncol(DLBCL1))
DLBCL2 <- matrix(ifelse(is.na(as.numeric(DLBCL2)), 0, as.numeric(DLBCL2)), nrow = nrow(DLBCL2), ncol = ncol(DLBCL2))
FL <- matrix(ifelse(is.na(as.numeric(FL)), 0, as.numeric(FL)), nrow = nrow(FL), ncol = ncol(FL))
CLL <- matrix(ifelse(is.na(as.numeric(CLL)), 0, as.numeric(CLL)), nrow = nrow(CLL), ncol = ncol(CLL))

# Output the dimensions of each matrix to verify
print(dim(DLBCL1))
print(dim(DLBCL2))
print(dim(FL))
print(dim(CLL))


group1 <- t(DLBCL1)
group2 <- t(DLBCL2)

dim(group1)
dim(group2)

library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# rats.xls
{
# Step 1: Load the necessary library to read the Excel file
library(readxl)

# Step 2: Read the excel file containing the rats data
rats_data <- read_excel("rats.xls")

# Step 3: Remove the GROUP and T_0 columns from the dataset
rats_data_clean <- rats_data[, -c(1, 2)]

# Step 4: Split the data into groups based on the GROUP column values (1, 2, 3, 4)
group1_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 1, ]))
group2_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 2, ]))
group3_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 3, ]))
group4_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 4, ]))

# Step 5: Output the dimensions of each group matrix to verify that the split is correct
print(dim(group1_matrix))
print(dim(group2_matrix))
print(dim(group3_matrix))
print(dim(group4_matrix))


group1 <- t(group1_matrix)
group2 <- t(group2_matrix)

dim(group1)
dim(group2)

library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# Pancreatic
{
# GSE24279_non-normalized.txt.gz
# Load the data from the compressed .gz file
data_GSE24279 <- read.delim(gzfile("GSE24279_non-normalized.txt.gz"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Inspect the dataset
# head(data_GSE24279)

# Create matrices for each class directly from the dataset
normal_controls <- as.matrix(data_GSE24279[, grep("healthy", colnames(data_GSE24279))])
pancreatitis <- as.matrix(data_GSE24279[, grep("pancreatitis", colnames(data_GSE24279))])
pancreatic_cancer <- as.matrix(data_GSE24279[, grep("pancreatic_cancer", colnames(data_GSE24279))])

# Output the dimensions of each matrix to verify
print(dim(normal_controls))
print(dim(pancreatitis))
print(dim(pancreatic_cancer))

pancreatic<- cbind(pancreatitis,pancreatic_cancer)

group1 <- t(pancreatic)
group2 <- t(normal_controls)

dim(group1)
dim(group2)

library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# Heart Disease
{
# processed_cleveland_data.csv
# Load the dataset
data_cleveland <- read.csv("processed_cleveland_data.csv", na.strings = "?")
data_cleveland[is.na(data_cleveland)] <- 0

# Transpose the dataset, setting target as the column names
data_cleveland_t <- t(data_cleveland[, 1:(ncol(data_cleveland)-1)])  # Transpose the features, excluding the target column
colnames(data_cleveland_t) <- data_cleveland$target  # Set target as the column names

# Verify the transposed dataset
# head(data_cleveland_t)

# Create matrices for each target class
class_0 <- as.matrix(data_cleveland_t[, grep("0", colnames(data_cleveland_t))])
class_1 <- as.matrix(data_cleveland_t[, grep("1", colnames(data_cleveland_t))])
class_2 <- as.matrix(data_cleveland_t[, grep("2", colnames(data_cleveland_t))])
class_3 <- as.matrix(data_cleveland_t[, grep("3", colnames(data_cleveland_t))])
class_4 <- as.matrix(data_cleveland_t[, grep("4", colnames(data_cleveland_t))])

is.numeric(class_0)
is.numeric(class_1)
is.numeric(class_2)
is.numeric(class_3)
is.numeric(class_4)

# Output the dimensions of each matrix to verify
print(dim(class_0))
print(dim(class_1))
print(dim(class_2))
print(dim(class_3))
print(dim(class_4))

group1 <- t(class_1)
group2 <- t(class_2)

dim(group1)
dim(group2)

library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## highDmean
# SKK_test(group1, group2)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT(group1, group2)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## SHT
# mean2.1996BS(group1,group2)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART(group1, group2)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD(group1,group2)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## highmean
# apval_Bai1996(group1,group2)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# apval_Sri2008(group1,group2)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT(group1,group2)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## Chen and Qin 2010
# apval_Chen2010(group1,group2, eq.cov = TRUE)
run_test_10_times_Chen2010 <- function(group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    apval_Chen2010(group1, group2, eq.cov = TRUE)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

average_time_apval_Chen2010 <- run_test_10_times_Chen2010(group1, group2)
cat("Average time for apval_Chen2010: ", average_time_apval_Chen2010, " seconds\n")

# CQ2010.TSBF.NABT(group1, group2)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}

########### time

# keletal muscle
{# GSE25941_series_matrix.txt
# Load the skeletal muscle dataset
# The dataset is in a series matrix format, and the metadata and expression matrix are extracted for analysis

# Step 1: Read the entire file
# Load the GSE25941_series_matrix.txt file into memory
file_content <- readLines("GSE25941_series_matrix.txt")

# Step 2: Extract metadata (Sample Characteristics)
# Extract sample characteristics for gender and age of each sample
sample_characteristics_ch2_line <- grep("!Sample_characteristics_ch2", file_content)
sample_characteristics_ch2 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch2_line]), "\t"))[-1]

sample_characteristics_ch3_line <- grep("!Sample_characteristics_ch3", file_content)
sample_characteristics_ch3 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch3_line]), "\t"))[-1]

# Step 3: Read the expression matrix
# Locate the start of the expression matrix and load the matrix data from the file
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
# Combine the metadata and the gene expression matrix for further analysis
combined_GSE25941 <- rbind(metadata_df, expression_matrix)

# Step 6: Create matrices for each group (young females, young males, old females, old males)
# Extract the data for each group based on gender and age information from the metadata
# Exclude metadata rows and the last row for consistent formatting
young_females_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Female" & combined_GSE25941[1, ] == "age: Young"])
young_males_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Male" & combined_GSE25941[1, ] == "age: Young"])
old_females_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Female" & combined_GSE25941[1, ] == "age: Old"])
old_males_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Male" & combined_GSE25941[1, ] == "age: Old"])

# Convert character data to numeric, handling potential non-numeric values
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
# Print the dimensions of each group matrix to ensure correct extraction and scaling
print(dim(young_females_matrix))
print(dim(young_males_matrix))
print(dim(old_females_matrix))
print(dim(old_males_matrix))

# Step 9: Combine young females and young males matrices, and old females and old males matrices
# Combine the matrices for young and old groups for further analysis
young_combined_matrix <- cbind(young_females_matrix, young_males_matrix)
old_combined_matrix <- cbind(old_females_matrix, old_males_matrix)

# Step 10: Output the dimensions of each combined matrix to verify correctness
# Print the dimensions of each combined group matrix to ensure correct combination
print(dim(young_combined_matrix))
print(dim(old_combined_matrix))

# Transpose the combined matrices to match the format for hypothesis testing
group1 <- t(young_combined_matrix)
group2 <- t(old_combined_matrix)

# Output the dimensions of the transposed matrices
dim(group1)
dim(group2)

# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# SARS-CoV-2
{# GSE156063_series_matrix.txt
# GSE156063_series_matrix.txt Processing Script

# Step 1: Read the series matrix file (metadata file)
# Load the GSE156063_series_matrix.txt file into memory
file_content <- readLines("GSE156063_series_matrix.txt")

# Step 2: Extract Sample Characteristics (Sample Characteristics ch2: disease states)
# Locate the lines containing disease states information for each sample
sample_characteristics_ch2_line <- grep("!Sample_characteristics_ch2", file_content)

# Split the characteristics data into individual elements by tab separator
# Remove any quotation marks and extract the characteristics data
sample_characteristics_ch2 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch2_line]), "\t"))[-1]

# Step 3: Read the gene expression data from the CSV file (GSE156063_swab_gene_counts.csv)
# Load the gene expression matrix from the CSV file
gene_expression_data <- read.csv("GSE156063_swab_gene_counts.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Step 4: Convert sample characteristics to rows
# Create a new data frame with sample characteristics ch2 as a row
# This data frame will store the metadata as the first row in the combined dataset
metadata_df <- data.frame(
  ID_REF = "Sample_Characteristics_Ch2",  # Label for the metadata row
  t(sample_characteristics_ch2),  # Transpose the sample characteristics into a row format
  stringsAsFactors = FALSE
)

# Step 5: Ensure that the column names of the metadata match the gene expression data
# Update the column names of the metadata to match the column names of the expression data for consistency
colnames(metadata_df) <- colnames(gene_expression_data)

# Step 6: Combine the metadata and the gene expression data
# Bind the metadata row to the gene expression matrix for further analysis
combined_GSE156063 <- rbind(metadata_df, gene_expression_data)

# Step 7: Save the combined data as a CSV file
# Save the combined metadata and expression data to a CSV file for record-keeping and subsequent use
# write.csv(combined_GSE156063, "GSE156063_combined_data.csv", row.names = FALSE)
# Output the dimensions of the combined data to verify the data structure
dim(combined_GSE156063)

# Step 8: Create matrices for each disease state group (SC2, other viral, no virus)
# Select samples corresponding to each group based on disease state information from metadata
sc2_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: SC2"])
other_virus_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: other virus"])
no_virus_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: no virus"])

# Step 9: Convert character data to numeric if necessary
# Convert the matrices from character to numeric to prepare for downstream analysis
sc2_matrix <- apply(sc2_matrix, 2, as.numeric)
other_virus_matrix <- apply(other_virus_matrix, 2, as.numeric)
no_virus_matrix <- apply(no_virus_matrix, 2, as.numeric)

# Step 10: Output the dimensions of each matrix to verify
# Print the dimensions of each matrix to ensure correct extraction of samples by disease state
print(dim(sc2_matrix))
print(dim(other_virus_matrix))
print(dim(no_virus_matrix))

# Step 11: Combine SC2 and other viral matrices into a new matrix
# Combine SC2 and other viral samples into a single matrix for further comparison with no virus group
sc2_other_combined_matrix <- cbind(sc2_matrix, other_virus_matrix)

# Step 12: Output the dimensions of the combined matrix and the no virus matrix to verify
# Print the dimensions of the combined SC2/other viral matrix and the no virus matrix to ensure correct combination
print(dim(sc2_other_combined_matrix))
print(dim(no_virus_matrix))

# Step 13: Transpose the combined matrices for downstream analysis
# Prepare the combined groups for hypothesis testing by transposing the matrices
group1 <- t(sc2_other_combined_matrix)
group2 <- t(no_virus_matrix)

# Output the dimensions of the transposed matrices to verify correctness
dim(group1)
dim(group2)


# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# COVID-19 from HDNRA
{
# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Step 2: Load the COVID-19 dataset
# Load the COVID-19 dataset from the HDNRA package, which contains gene expression data for COVID-19 patients
data("COVID19")

# Step 3: Extract group data for hypothesis testing
# Extract specific rows from the dataset to create two groups for comparison

# Create 'group1' containing specific rows (samples 2 to 19 and 82 to 87) representing certain subsets of the data
group1 <- as.matrix(COVID19[c(2:19, 82:87), ])
# Output the dimensions of group1 to verify the extraction
dim(group1)

# Create 'group2' containing the remaining rows (excluding samples 1 to 19 and 82 to 87)
group2 <- as.matrix(COVID19[-c(1:19, 82:87), ])
# Output the dimensions of group2 to verify the extraction
dim(group2)

# Step 4: Verify the dimensions of each group
# Print the dimensions of both groups to ensure that the data was correctly divided into group1 and group2
dim(group1)
dim(group2)


# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# yeoh-2002-v2
{
# Step 1: Load the Yeoh 2002 dataset
# Load the dataset from the text file and skip the first row
# This dataset contains gene expression information for different leukemia subtypes
data1 <- read.table("yeoh-2002-v2_database.txt", header = TRUE, sep = "\t", fill = TRUE, skip = 1)

# Step 2: Verify the data structure
# Uncomment the lines below to examine the first few rows and column names of the dataset
# This is helpful for understanding the data format and verifying the correct loading of data
# head(data1)
# colnames(data1)

# Step 3: Create matrices for each sample type
# Extract samples of different leukemia subtypes into separate matrices for analysis

# Extract columns corresponding to 'BCR-ABL' samples
BCR_ABL <- as.matrix(data1[, grep("BCR.ABL", names(data1))])

# Extract columns corresponding to 'E2A-PBX1' samples
E2A_PBX1 <- as.matrix(data1[, grep("E2A.PBX1", names(data1))])

# Extract columns corresponding to 'Hyperdiploid 50' samples
Hyperdiploid_50 <- as.matrix(data1[, grep("Hyperdip.50", names(data1))])

# Extract columns corresponding to 'MLL' samples
MLL <- as.matrix(data1[, grep("MLL", names(data1))])

# Extract columns corresponding to 'TEL-AML1' samples
TEL_AML1 <- as.matrix(data1[, grep("TEL.AML", names(data1))])

# Extract columns corresponding to 'T-ALL' samples
T_ALL <- as.matrix(data1[, grep("T.ALL", names(data1))])

# Step 4: Output the dimensions of each matrix to verify
# Print the dimensions of each matrix to ensure correct extraction of sample data
print(dim(BCR_ABL))
print(dim(E2A_PBX1))
print(dim(Hyperdiploid_50))
print(dim(MLL))
print(dim(TEL_AML1))
print(dim(T_ALL))

# Step 5: Create group matrices for downstream analysis
# Transpose matrices for selected groups to prepare for further hypothesis testing or comparison
# Create 'group1' from the 'BCR-ABL' samples
group1 <- t(BCR_ABL)

# Create 'group2' from the 'TEL-AML1' samples
group2 <- t(TEL_AML1)

# Step 6: Output the dimensions of the transposed matrices to verify
# Print the dimensions of each group matrix to ensure they were correctly transposed
dim(group1)
dim(group2)

# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# alizadeh-2000-v3
{
# Step 1: Load the Alizadeh 2000 dataset
# Load the dataset from the text file and skip the first row
# This dataset contains gene expression information for different lymphoma subtypes
data2 <- read.table("alizadeh-2000-v3_database.txt", header = TRUE, sep = "\t", fill = TRUE, skip = 1)

# Step 2: Verify the data structure
# Uncomment the line below to examine the first few rows of the dataset
# This is helpful for understanding the data format and verifying the correct loading of data
# head(data2)

# Step 3: Create matrices for each sample type
# Extract samples of different lymphoma subtypes into separate matrices for analysis

# Extract columns corresponding to 'DLBCL1' samples
DLBCL1 <- as.matrix(data2[, grep("DLBCL1", names(data2))])

# Extract columns corresponding to 'DLBCL2' samples
DLBCL2 <- as.matrix(data2[, grep("DLBCL2", names(data2))])

# Extract columns corresponding to 'FL' (Follicular Lymphoma) samples
FL <- as.matrix(data2[, grep("FL", names(data2))])

# Extract columns corresponding to 'CLL' (Chronic Lymphocytic Leukemia) samples
CLL <- as.matrix(data2[, grep("CLL", names(data2))])

# Step 4: Convert the matrices to numeric, handling any non-numeric values
# Convert all data to numeric values for further analysis

# Convert the 'DLBCL1' matrix to numeric, replacing non-numeric values with NA
DLBCL1 <- apply(DLBCL1, 2, as.numeric)

# Convert the 'DLBCL2' matrix to numeric, replacing non-numeric values with NA
DLBCL2 <- apply(DLBCL2, 2, as.numeric)

# Convert the 'FL' matrix to numeric, replacing non-numeric values with NA
FL <- apply(FL, 2, as.numeric)

# Convert the 'CLL' matrix to numeric, replacing non-numeric values with NA
CLL <- apply(CLL, 2, as.numeric)

# Step 5: Handle NA values by replacing them with 0
# Convert the numeric matrices to handle NA values by replacing them with 0

# Replace NA values in the 'DLBCL1' matrix with 0
DLBCL1 <- matrix(ifelse(is.na(as.numeric(DLBCL1)), 0, as.numeric(DLBCL1)), nrow = nrow(DLBCL1), ncol = ncol(DLBCL1))

# Replace NA values in the 'DLBCL2' matrix with 0
DLBCL2 <- matrix(ifelse(is.na(as.numeric(DLBCL2)), 0, as.numeric(DLBCL2)), nrow = nrow(DLBCL2), ncol = ncol(DLBCL2))

# Replace NA values in the 'FL' matrix with 0
FL <- matrix(ifelse(is.na(as.numeric(FL)), 0, as.numeric(FL)), nrow = nrow(FL), ncol = ncol(FL))

# Replace NA values in the 'CLL' matrix with 0
CLL <- matrix(ifelse(is.na(as.numeric(CLL)), 0, as.numeric(CLL)), nrow = nrow(CLL), ncol = ncol(CLL))

# Step 6: Output the dimensions of each matrix to verify
# Print the dimensions of each matrix to ensure correct extraction and conversion
print(dim(DLBCL1))
print(dim(DLBCL2))
print(dim(FL))
print(dim(CLL))

# Step 7: Create group matrices for downstream analysis
# Transpose matrices for selected groups to prepare for further hypothesis testing or comparison

# Create 'group1' from the 'DLBCL1' samples
group1 <- t(DLBCL1)

# Create 'group2' from the 'DLBCL2' samples
group2 <- t(DLBCL2)

# Step 8: Output the dimensions of the transposed matrices to verify
# Print the dimensions of each group matrix to ensure they were correctly transposed
dim(group1)
dim(group2)


# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# rats.xls
{
# Step 1: Load the necessary library to read the Excel file
# Load the 'readxl' library to read data from an Excel file into R
library(readxl)

# Step 2: Read the Excel file containing the rats data
# Load the rats dataset from the 'rats.xls' Excel file
# This dataset contains information on multiple experimental groups of rats
rats_data <- read_excel("rats.xls")

# Step 3: Remove the GROUP and T_0 columns from the dataset
# Remove the first two columns (GROUP and T_0) from the dataset for analysis
# This leaves only the columns that contain time series data for further analysis
rats_data_clean <- rats_data[, -c(1, 2)]

# Step 4: Split the data into groups based on the GROUP column values (1, 2, 3, 4)
# Create matrices for each experimental group based on the GROUP column
# Extract rows corresponding to each group and convert them to matrices

# Extract and transpose data for 'Group 1' rats
group1_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 1, ]))

# Extract and transpose data for 'Group 2' rats
group2_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 2, ]))

# Extract and transpose data for 'Group 3' rats
group3_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 3, ]))

# Extract and transpose data for 'Group 4' rats
group4_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 4, ]))

# Step 5: Output the dimensions of each group matrix to verify that the split is correct
# Print the dimensions of each group matrix to ensure the correct assignment of samples to groups
print(dim(group1_matrix))
print(dim(group2_matrix))
print(dim(group3_matrix))
print(dim(group4_matrix))

# Step 6: Prepare group matrices for downstream analysis
# Transpose the matrices to prepare for further hypothesis testing or comparison
# Here, we take group1 and group2 for further analysis

group1 <- t(group1_matrix)
group2 <- t(group2_matrix)

# Step 7: Output the dimensions of the transposed matrices to verify correctness
# Print the dimensions of group1 and group2 to ensure they were correctly transposed
dim(group1)
dim(group2)


# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# Pancreatic
{ # GSE24279_non-normalized.txt.gz
# Step 1: Load the data from the compressed .gz file
# Load the dataset from the compressed text file (GSE24279_non-normalized.txt.gz)
# This dataset contains gene expression data for different classes of samples (healthy, pancreatitis, pancreatic cancer)
data_GSE24279 <- read.delim(gzfile("GSE24279_non-normalized.txt.gz"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 2: Inspect the dataset
# Uncomment the line below to examine the first few rows of the dataset
# This is useful for understanding the data structure and verifying the correct loading of data
# head(data_GSE24279)

# Step 3: Create matrices for each class directly from the dataset
# Extract samples of different classes into separate matrices for analysis

# Extract columns corresponding to 'healthy' samples (normal controls)
normal_controls <- as.matrix(data_GSE24279[, grep("healthy", colnames(data_GSE24279))])

# Extract columns corresponding to 'pancreatitis' samples
pancreatitis <- as.matrix(data_GSE24279[, grep("pancreatitis", colnames(data_GSE24279))])

# Extract columns corresponding to 'pancreatic cancer' samples
pancreatic_cancer <- as.matrix(data_GSE24279[, grep("pancreatic_cancer", colnames(data_GSE24279))])

# Step 4: Output the dimensions of each matrix to verify
# Print the dimensions of each matrix to ensure correct extraction of sample data by class
print(dim(normal_controls))
print(dim(pancreatitis))
print(dim(pancreatic_cancer))

# Step 5: Combine matrices for pancreatitis and pancreatic cancer
# Create a combined matrix for the 'pancreatitis' and 'pancreatic cancer' samples for further analysis
pancreatic <- cbind(pancreatitis, pancreatic_cancer)

# Step 6: Prepare group matrices for downstream analysis
# Transpose the matrices for selected groups (pancreatic samples and healthy controls) to prepare for hypothesis testing or comparison
group1 <- t(pancreatic)
group2 <- t(normal_controls)

# Step 7: Output the dimensions of the transposed matrices to verify correctness
# Print the dimensions of group1 (pancreatic samples) and group2 (healthy controls) to ensure they were correctly transposed
dim(group1)
dim(group2)

# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}


# Heart Disease
{# processed_cleveland_data.csv
# Step 1: Load the Cleveland heart disease dataset
# Load the dataset from the CSV file and replace missing values represented by '?' with 0
# This dataset contains medical records used for predicting heart disease
data_cleveland <- read.csv("processed_cleveland_data.csv", na.strings = "?")

# Replace NA values with 0 for further analysis
data_cleveland[is.na(data_cleveland)] <- 0

# Step 2: Transpose the dataset, setting target as the column names
# Transpose the features of the dataset, excluding the target column
# Set the target column values as the column names in the transposed dataset
data_cleveland_t <- t(data_cleveland[, 1:(ncol(data_cleveland) - 1)])  # Transpose the features
colnames(data_cleveland_t) <- data_cleveland$target  # Set target values as column names

# Step 3: Verify the transposed dataset
# Uncomment the line below to examine the first few rows of the transposed dataset
# This helps verify that the transposition and column assignments were done correctly
# head(data_cleveland_t)

# Step 4: Create matrices for each target class
# Extract samples corresponding to each target class into separate matrices

# Extract columns corresponding to 'class 0' (no heart disease)
class_0 <- as.matrix(data_cleveland_t[, grep("0", colnames(data_cleveland_t))])

# Extract columns corresponding to 'class 1' (mild heart disease)
class_1 <- as.matrix(data_cleveland_t[, grep("1", colnames(data_cleveland_t))])

# Extract columns corresponding to 'class 2' (moderate heart disease)
class_2 <- as.matrix(data_cleveland_t[, grep("2", colnames(data_cleveland_t))])

# Extract columns corresponding to 'class 3' (severe heart disease)
class_3 <- as.matrix(data_cleveland_t[, grep("3", colnames(data_cleveland_t))])

# Extract columns corresponding to 'class 4' (very severe heart disease)
class_4 <- as.matrix(data_cleveland_t[, grep("4", colnames(data_cleveland_t))])

# Step 5: Verify if the matrices contain numeric data
# Ensure that each class matrix contains only numeric data for downstream analysis
is.numeric(class_0)
is.numeric(class_1)
is.numeric(class_2)
is.numeric(class_3)
is.numeric(class_4)

# Step 6: Output the dimensions of each matrix to verify
# Print the dimensions of each class matrix to ensure correct extraction
print(dim(class_0))
print(dim(class_1))
print(dim(class_2))
print(dim(class_3))
print(dim(class_4))

# Step 7: Prepare group matrices for downstream analysis
# Transpose the matrices for selected groups to prepare for hypothesis testing or comparison
# Here, group1 corresponds to 'class 1' and group2 corresponds to 'class 2'

group1 <- t(class_1)
group2 <- t(class_2)

# Step 8: Output the dimensions of the transposed matrices to verify correctness
# Print the dimensions of group1 and group2 to ensure they were correctly transposed
dim(group1)
dim(group2)


# Load the necessary libraries for hypothesis testing
library(HDNRA)
library(highmean)
library(highDmean)
library(SHT)

# Function to run and time a test 10 times
# Helper function to compute the average time taken by each test over 10 repetitions
run_test_10_times <- function(test_function, group1, group2) {
  times <- replicate(10, {
    startTime <- Sys.time()
    test_function(group1, group2)
    endTime <- Sys.time()
    as.numeric(difftime(endTime, startTime, units = "secs"))
  })
  mean(times)  # Return the average time
}

## HDNRA and highDmean comparison
# SKK_test (high-dimensional two-sample test for mean differences)
average_time_SKK_test <- run_test_10_times(SKK_test, group1, group2)
cat("Average time for SKK_test: ", average_time_SKK_test, " seconds\n")

# SKK2013.TSBF.NABT (high-dimensional two-sample test based on Bai and Saranadasa's method)
average_time_SKK2013 <- run_test_10_times(SKK2013.TSBF.NABT, group1, group2)
cat("Average time for SKK2013.TSBF.NABT: ", average_time_SKK2013, " seconds\n")

## HDNRA and SHT comparison
# mean2.1996BS (Bai and Saranadasa's mean test)
average_time_mean2_1996BS <- run_test_10_times(mean2.1996BS, group1, group2)
cat("Average time for mean2.1996BS: ", average_time_mean2_1996BS, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# mean2.2008SD (Srivastava and Du's mean test)
average_time_mean2_2008SD <- run_test_10_times(mean2.2008SD, group1, group2)
cat("Average time for mean2.2008SD: ", average_time_mean2_2008SD, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

## HDNRA and highmean comparison
# apval_Bai1996 (Adjusted p-value using Bai's method)
average_time_apval_Bai1996 <- run_test_10_times(apval_Bai1996, group1, group2)
cat("Average time for apval_Bai1996: ", average_time_apval_Bai1996, " seconds\n")

# BS1996.TS.NART (Test statistic based on Bai and Saranadasa)
average_time_BS1996_TS_NART <- run_test_10_times(BS1996.TS.NART, group1, group2)
cat("Average time for BS1996.TS.NART: ", average_time_BS1996_TS_NART, " seconds\n")

# apval_Sri2008 (Adjusted p-value using Srivastava's method)
average_time_apval_Sri2008 <- run_test_10_times(apval_Sri2008, group1, group2)
cat("Average time for apval_Sri2008: ", average_time_apval_Sri2008, " seconds\n")

# SD2008.TS.NABT (Srivastava and Du's test statistic)
average_time_SD2008_TS_NABT <- run_test_10_times(SD2008.TS.NABT, group1, group2)
cat("Average time for SD2008.TS.NABT: ", average_time_SD2008_TS_NABT, " seconds\n")

# apval_Chen2010 (Chen and Qin's mean test for equal covariance)
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

# CQ2010.TSBF.NABT (Test statistic based on Chen and Qin's method)
average_time_CQ2010_TSBF_NABT <- run_test_10_times(CQ2010.TSBF.NABT, group1, group2)
cat("Average time for CQ2010.TSBF.NABT: ", average_time_CQ2010_TSBF_NABT, " seconds\n")
}

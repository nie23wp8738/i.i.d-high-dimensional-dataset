### high-dimensional datasets used by HDNRA package

## comments for the functions in HDNRA
# Perform hypothesis tests using normal-reference approach-based methods
# These functions perform GLHT using different test statistics and reference distributions
# Normal-reference approach based tests for GLHT problem
# ZGZ2017.GLHT.2cNRT(Y, G, n, p)  # Normal-reference approach based GLHT test using ZGZ 2017 method
# ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)  # Normal-reference approach based GLHT test using ZZZ 2022 method
# ZZG2022.GLHTBF.2cNRT(Y, G, n, p)  # Bayesian factor-based GLHT test using ZZG 2022 method
# ZZ2022.GLHTBF.3cNRT(Y, G, n, p)  # Bayesian factor-based GLHT test for 3-way comparison using ZZ 2022 method
# ZZ2022.GLHT.3cNRT(Y, G, n, p)  # Normal-reference approach based 3-way GLHT test using ZZ 2022 method

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
# FHW2004.GLHT.NABT(Y, X, C, n, p)  # Non-asymptotic bootstrap test for GLHT using FHW 2004 method
# SF2006.GLHT.NABT(Y, X, C, n, p)  # Non-asymptotic bootstrap test for GLHT using SF 2006 method
# S2007.ks.NABT(Y, n, p)  # Kolmogorov-Smirnov test for GLHT using S 2007 method
# YS2012.GLHT.NABT(Y, X, C, n, p)  # Non-asymptotic bootstrap test for GLHT using YS 2012 method
# ZGZ2017.GLHTBF.NABT(Y, G, n, p)  # Bayesian factor-based non-asymptotic bootstrap test for GLHT using ZGZ 2017 method


# data 1
{# Load the Yeoh 2002 dataset and skip the first row as it contains metadata
# The dataset is tab-separated and the header is present
# data1 contains gene expression data for multiple leukemia subtypes
data1 <- read.table("yeoh-2002-v2_database.txt", header=TRUE, sep="\t", fill=TRUE, skip=1)

# Verify the data structure by examining the first few rows and column names
# head(data1)
# colnames(data1)

# Create matrices for each leukemia subtype using pattern matching on column names
BCR_ABL <- as.matrix(data1[, grep("BCR.ABL", names(data1))])
E2A_PBX1 <- as.matrix(data1[, grep("E2A.PBX1", names(data1))])
Hyperdiploid_50 <- as.matrix(data1[, grep("Hyperdip.50", names(data1))])
MLL <- as.matrix(data1[, grep("MLL", names(data1))])
TEL_AML1 <- as.matrix(data1[, grep("TEL.AML", names(data1))])
T_ALL <- as.matrix(data1[, grep("T.ALL", names(data1))])

# Output the dimensions of each matrix to verify that data has been loaded correctly
print(dim(BCR_ABL))
print(dim(E2A_PBX1))
print(dim(Hyperdiploid_50))
print(dim(MLL))
print(dim(TEL_AML1))
print(dim(T_ALL))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of genes (features) from the T_ALL matrix
p <- dim(T_ALL)[1]

# Create a list of matrices, where each element is a transposed matrix of gene expression data
# for each leukemia subtype
Y <- list()
k <- 6
Y[[1]] <- t(BCR_ABL)
Y[[2]] <- t(E2A_PBX1)
Y[[3]] <- t(Hyperdiploid_50)
Y[[4]] <- t(MLL)
Y[[5]] <- t(TEL_AML1)
Y[[6]] <- t(T_ALL)

# Define the number of samples for each leukemia subtype
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]), nrow(Y[[5]]), nrow(Y[[6]]))

# Create the contrast matrix for group comparisons (General Linear Hypothesis Testing setup)
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for GLHT problem
# Define the dimension of hypothesis test parameters
q <- k - 1

# Construct the design matrix for hypothesis testing (involving comparisons between groups)
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), rep(1, n[3]), rep(0, sum(n)), rep(1, n[4]), rep(0, sum(n)), rep(1, n[5]), rep(0, sum(n)), rep(1, n[6])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach based methods
# Normal-reference approach based test for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# Other tests for GLHT problem
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 2
{# Load the Alizadeh 2000 dataset
# Skip the first row as it contains metadata
# The dataset is tab-separated and the header is present
data2 <- read.table("alizadeh-2000-v3_database.txt", header=TRUE, sep="\t", fill=TRUE, skip=1)

# Verify the data structure by examining the first few rows
# head(data2)

# Create matrices for each sample type using pattern matching on column names
DLBCL1 <- as.matrix(data2[, grep("DLBCL1", names(data2))])
DLBCL2 <- as.matrix(data2[, grep("DLBCL2", names(data2))])
FL <- as.matrix(data2[, grep("FL", names(data2))])
CLL <- as.matrix(data2[, grep("CLL", names(data2))])

# Convert the matrices to numeric values, handling any non-numeric values
# This step ensures all values are numeric for further analysis
DLBCL1 <- apply(DLBCL1, 2, as.numeric)
DLBCL2 <- apply(DLBCL2, 2, as.numeric)
FL <- apply(FL, 2, as.numeric)
CLL <- apply(CLL, 2, as.numeric)

# Replace any NA values with 0 to avoid issues in analysis
DLBCL1 <- matrix(ifelse(is.na(as.numeric(DLBCL1)), 0, as.numeric(DLBCL1)), nrow = nrow(DLBCL1), ncol = ncol(DLBCL1))
DLBCL2 <- matrix(ifelse(is.na(as.numeric(DLBCL2)), 0, as.numeric(DLBCL2)), nrow = nrow(DLBCL2), ncol = ncol(DLBCL2))
FL <- matrix(ifelse(is.na(as.numeric(FL)), 0, as.numeric(FL)), nrow = nrow(FL), ncol = ncol(FL))
CLL <- matrix(ifelse(is.na(as.numeric(CLL)), 0, as.numeric(CLL)), nrow = nrow(CLL), ncol = ncol(CLL))

# Output the dimensions of each matrix to verify that data has been loaded correctly
print(dim(DLBCL1))
print(dim(DLBCL2))
print(dim(FL))
print(dim(CLL))

# Check if there are any remaining NA values
# This step is for verification purposes to ensure data integrity
# cat("DLBCL1 NA count:", sum(is.na(DLBCL1)), "\n")
# cat("DLBCL2 NA count:", sum(is.na(DLBCL2)), "\n")
# cat("FL NA count:", sum(is.na(FL)), "\n")
# cat("CLL NA count:", sum(is.na(CLL)), "\n")
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of genes (features) from the DLBCL1 matrix
p <- dim(DLBCL1)[1]

# Create a list of matrices, where each element is a transposed matrix of gene expression data for each sample type
Y <- list()
k <- 4
Y[[1]] <- t(DLBCL1)
Y[[2]] <- t(DLBCL2)
Y[[3]] <- t(FL)
Y[[4]] <- t(CLL)

# Define the number of samples for each sample type
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]))

# Create the contrast matrix for group comparisons (General Linear Hypothesis Testing setup)
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for GLHT problem
# Define the dimension of hypothesis test parameters
q <- k - 1

# Construct the design matrix for hypothesis testing (involving comparisons between groups)
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), rep(1, n[3]), rep(0, sum(n)), rep(1, n[4])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach based methods
# Normal-reference approach based test for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# Other tests for GLHT problem
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 3
{
# Load the Tomlins 2006 dataset, skipping the first row containing metadata
# The dataset is tab-separated and includes column headers
data3 <- read.table("tomlins-2006_database.txt", header=TRUE, sep="\t", fill=TRUE, skip=1)

# Verify the data structure by examining the first few rows and column names
# Uncomment the lines below to check the structure of data3
# head(data3)  # View the first few rows of the dataset
# colnames(data3)  # View the column names

# Extract matrices for each sample type using pattern matching on column names
# EPI, MET, PCA, PIN, and STROMA represent different prostate tissue types
EPI <- as.matrix(data3[, grep("EPI", names(data3))])
MET <- as.matrix(data3[, grep("MET", names(data3))])
PCA <- as.matrix(data3[, grep("PCA", names(data3))])
PIN <- as.matrix(data3[, grep("PIN", names(data3))])
STROMA <- as.matrix(data3[, grep("STROMA", names(data3))])

# Output the dimensions of each matrix to verify correct data extraction
print(dim(EPI))
print(dim(MET))
print(dim(PCA))
print(dim(PIN))
print(dim(STROMA))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of genes (features) based on the first matrix (EPI)
# Assuming all matrices have the same number of genes, use the number of rows from EPI
p <- dim(EPI)[1]

# Create a list of matrices, each containing the transposed gene expression data for each sample type
# Transposing to match the expected format for hypothesis testing (genes as columns)
Y <- list()
k <- 5  # Total number of sample types
Y[[1]] <- t(EPI)
Y[[2]] <- t(MET)
Y[[3]] <- t(PCA)
Y[[4]] <- t(PIN)
Y[[5]] <- t(STROMA)

# Define the number of samples for each sample type (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]), nrow(Y[[5]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4]), rep(0, sum(n)), 
              rep(1, n[5])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)

# data 4
{
# Load the Iris dataset
# The dataset is in CSV format and includes column headers
data_iris <- read.csv("Iris Data.csv")

# Transpose the dataset, keeping only the features for transposition
# Set the species names as the column names after transposition
data_iris_t <- t(data_iris[, 1:4])  # Transpose the feature columns
colnames(data_iris_t) <- data_iris$species  # Set species names as column names

# Verify the transposed dataset by examining the first few rows
# Uncomment the line below to check the transposed dataset
# head(data_iris_t)

# Extract matrices for each species using pattern matching on column names
# Setosa, Versicolor, and Virginica represent the three species in the Iris dataset
Setosa <- as.matrix(data_iris_t[, grep("setosa", colnames(data_iris_t), ignore.case = TRUE)])
Versicolor <- as.matrix(data_iris_t[, grep("versicolor", colnames(data_iris_t), ignore.case = TRUE)])
Virginica <- as.matrix(data_iris_t[, grep("virginica", colnames(data_iris_t), ignore.case = TRUE)])

# Output the dimensions of each matrix to verify correct data extraction
print(dim(Setosa))
print(dim(Versicolor))
print(dim(Virginica))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (Setosa)
# Assuming all matrices have the same number of features, use the number of rows from Setosa
p <- dim(Setosa)[1]

# Create a list of matrices, each containing the transposed feature data for each species
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 3  # Total number of species
Y[[1]] <- t(Setosa)
Y[[2]] <- t(Versicolor)
Y[[3]] <- t(Virginica)

# Define the number of samples for each species (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 5
{# ecoli_data.csv
# Load the E. coli dataset
# The dataset is in CSV format and includes column headers
data_ecoli <- read.csv("ecoli_data.csv")

# Transpose the dataset, keeping only the features for transposition
# Set the class labels as the column names after transposition
data_ecoli_t <- t(data_ecoli[, 1:7])  # Transpose the feature columns
colnames(data_ecoli_t) <- data_ecoli$class  # Set class labels as column names

# Verify the transposed dataset by examining the first few rows
# Uncomment the line below to check the transposed dataset
# head(data_ecoli_t)

# Extract matrices for each class using pattern matching on column names
# cp, im, pp, imU, om, omL, imS, and imL represent different localization sites in the E. coli dataset
cp <- as.matrix(data_ecoli_t[, grep("cp", colnames(data_ecoli_t), ignore.case = TRUE)])
im <- as.matrix(data_ecoli_t[, grep("im", colnames(data_ecoli_t), ignore.case = TRUE)])  # Match "im" but not "imL", "imU", "imS"
pp <- as.matrix(data_ecoli_t[, grep("pp", colnames(data_ecoli_t), ignore.case = TRUE)])
imU <- as.matrix(data_ecoli_t[, grep("imU", colnames(data_ecoli_t), ignore.case = TRUE)])
om <- as.matrix(data_ecoli_t[, grep("om", colnames(data_ecoli_t), ignore.case = TRUE)])  # Match "om" but not "omL"
omL <- as.matrix(data_ecoli_t[, grep("omL", colnames(data_ecoli_t), ignore.case = TRUE)])
imS <- as.matrix(data_ecoli_t[, grep("imS", colnames(data_ecoli_t), ignore.case = TRUE)])
imL <- as.matrix(data_ecoli_t[, grep("imL", colnames(data_ecoli_t), ignore.case = TRUE)])

# Output the dimensions of each matrix to verify correct data extraction
print(dim(cp))
print(dim(im))
print(dim(pp))
print(dim(imU))
print(dim(om))
print(dim(omL))
print(dim(imS))
print(dim(imL))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (cp)
# Assuming all matrices have the same number of features, use the number of rows from cp
p <- dim(cp)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 5  # Total number of classes
Y[[1]] <- t(cp)
Y[[2]] <- t(im)
Y[[3]] <- t(pp)
Y[[4]] <- t(imU)
Y[[5]] <- t(om)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]), nrow(Y[[5]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4]), rep(0, sum(n)), 
              rep(1, n[5])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 6
{# yeast_data.csv
# Load the Yeast dataset
# The dataset is in CSV format and includes column headers
data_yeast <- read.csv("yeast_data.csv")

# Transpose the dataset, keeping only the features for transposition
# Set the class labels as the column names after transposition
data_yeast_t <- t(data_yeast[, 1:8])  # Transpose the feature columns
colnames(data_yeast_t) <- data_yeast$class  # Set class labels as column names

# Verify the transposed dataset by examining the first few rows
# Uncomment the line below to check the transposed dataset
# head(data_yeast_t)

# Extract matrices for each class using pattern matching on column names
# CYT, NUC, MIT, ME3, ME2, ME1, EXC, VAC, POX, ERL represent different localization sites in the yeast dataset
CYT <- as.matrix(data_yeast_t[, grep("CYT", colnames(data_yeast_t), ignore.case = TRUE)])
NUC <- as.matrix(data_yeast_t[, grep("NUC", colnames(data_yeast_t), ignore.case = TRUE)])
MIT <- as.matrix(data_yeast_t[, grep("MIT", colnames(data_yeast_t), ignore.case = TRUE)])
ME3 <- as.matrix(data_yeast_t[, grep("ME3", colnames(data_yeast_t), ignore.case = TRUE)])
ME2 <- as.matrix(data_yeast_t[, grep("ME2", colnames(data_yeast_t), ignore.case = TRUE)])
ME1 <- as.matrix(data_yeast_t[, grep("ME1", colnames(data_yeast_t), ignore.case = TRUE)])
EXC <- as.matrix(data_yeast_t[, grep("EXC", colnames(data_yeast_t), ignore.case = TRUE)])
VAC <- as.matrix(data_yeast_t[, grep("VAC", colnames(data_yeast_t), ignore.case = TRUE)])
POX <- as.matrix(data_yeast_t[, grep("POX", colnames(data_yeast_t), ignore.case = TRUE)])
ERL <- as.matrix(data_yeast_t[, grep("ERL", colnames(data_yeast_t), ignore.case = TRUE)])

# Output the dimensions of each matrix to verify correct data extraction
print(dim(CYT))
print(dim(NUC))
print(dim(MIT))
print(dim(ME3))
print(dim(ME2))
print(dim(ME1))
print(dim(EXC))
print(dim(VAC))
print(dim(POX))
print(dim(ERL))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (CYT)
# Assuming all matrices have the same number of features, use the number of rows from CYT
p <- dim(CYT)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 9  # Total number of classes
Y[[1]] <- t(CYT)
Y[[2]] <- t(NUC)
Y[[3]] <- t(MIT)
Y[[4]] <- t(ME3)
Y[[5]] <- t(ME2)
Y[[6]] <- t(ME1)
Y[[7]] <- t(EXC)
Y[[8]] <- t(VAC)
Y[[9]] <- t(POX)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]), nrow(Y[[5]]), nrow(Y[[6]]), nrow(Y[[7]]), nrow(Y[[8]]), nrow(Y[[9]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4]), rep(0, sum(n)), 
              rep(1, n[5]), rep(0, sum(n)), rep(1, n[6]), rep(0, sum(n)), 
              rep(1, n[7]), rep(0, sum(n)), rep(1, n[8]), rep(0, sum(n)), 
              rep(1, n[9])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 7 Heart Disease
{# processed_cleveland_data.csv
# Load the Heart Disease dataset (Cleveland data)
# The dataset is in CSV format, and missing values are replaced with 0
data_cleveland <- read.csv("processed_cleveland_data.csv", na.strings = "?")
data_cleveland[is.na(data_cleveland)] <- 0  # Replace missing values with 0

# Transpose the dataset, keeping only the features for transposition
# Set the target labels as the column names after transposition
data_cleveland_t <- t(data_cleveland[, 1:(ncol(data_cleveland) - 1)])  # Transpose the feature columns, excluding the target column
colnames(data_cleveland_t) <- data_cleveland$target  # Set target labels as column names

# Verify the transposed dataset by examining the first few rows
# Uncomment the line below to check the transposed dataset
# head(data_cleveland_t)

# Extract matrices for each target class using pattern matching on column names
# Class labels 0, 1, 2, 3, 4 represent different levels of heart disease severity
class_0 <- as.matrix(data_cleveland_t[, grep("0", colnames(data_cleveland_t))])
class_1 <- as.matrix(data_cleveland_t[, grep("1", colnames(data_cleveland_t))])
class_2 <- as.matrix(data_cleveland_t[, grep("2", colnames(data_cleveland_t))])
class_3 <- as.matrix(data_cleveland_t[, grep("3", colnames(data_cleveland_t))])
class_4 <- as.matrix(data_cleveland_t[, grep("4", colnames(data_cleveland_t))])

# Verify that all extracted matrices are numeric
is.numeric(class_0)
is.numeric(class_1)
is.numeric(class_2)
is.numeric(class_3)
is.numeric(class_4)

# Output the dimensions of each matrix to verify correct data extraction
print(dim(class_0))
print(dim(class_1))
print(dim(class_2))
print(dim(class_3))
print(dim(class_4))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (class_0)
# Assuming all matrices have the same number of features, use the number of rows from class_0
p <- dim(class_0)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 5  # Total number of classes
Y[[1]] <- t(class_0)
Y[[2]] <- t(class_1)
Y[[3]] <- t(class_2)
Y[[4]] <- t(class_3)
Y[[5]] <- t(class_4)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]), nrow(Y[[5]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4]), rep(0, sum(n)), 
              rep(1, n[5])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 8
{# wine.csv
# Load the Wine dataset
# The dataset is in CSV format and includes class labels in the first column
data_wine <- read.csv("wine.csv")

# Transpose the dataset, keeping only the features for transposition
# Set the class labels as the column names after transposition
data_wine_t <- t(data_wine[, 2:ncol(data_wine)])  # Transpose the feature columns, excluding the class column
colnames(data_wine_t) <- data_wine[, 1]  # Set the class column as the column names

# Verify the transposed dataset by examining the first few rows
# Uncomment the line below to check the transposed dataset
# head(data_wine_t)

# Extract matrices for each class using pattern matching on column names
# Class labels 1, 2, 3 represent different types of wine
class_11 <- as.matrix(data_wine_t[, grep("1", colnames(data_wine_t))])
class_22 <- as.matrix(data_wine_t[, grep("2", colnames(data_wine_t))])
class_33 <- as.matrix(data_wine_t[, grep("3", colnames(data_wine_t))])

# Output the dimensions of each matrix to verify correct data extraction
print(dim(class_11))
print(dim(class_22))
print(dim(class_33))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (class_11)
# Assuming all matrices have the same number of features, use the number of rows from class_11
p <- dim(class_11)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 3  # Total number of classes
Y[[1]] <- t(class_11)
Y[[2]] <- t(class_22)
Y[[3]] <- t(class_33)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 9 Pancreatic
{# Load the Pancreatic dataset (GSE24279)
# The dataset is in a compressed .gz file and is tab-separated
# Load the data from the compressed file
# Missing values (if any) are not explicitly handled here
data_GSE24279 <- read.delim(gzfile("GSE24279_non-normalized.txt.gz"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Inspect the dataset by examining the first few rows
# Uncomment the line below to check the initial few rows
# head(data_GSE24279)

# Extract matrices for each class using pattern matching on column names
# Classes are 'healthy', 'pancreatitis', and 'pancreatic_cancer'
normal_controls <- as.matrix(data_GSE24279[, grep("healthy", colnames(data_GSE24279))])
pancreatitis <- as.matrix(data_GSE24279[, grep("pancreatitis", colnames(data_GSE24279))])
pancreatic_cancer <- as.matrix(data_GSE24279[, grep("pancreatic_cancer", colnames(data_GSE24279))])

# Output the dimensions of each matrix to verify correct data extraction
print(dim(normal_controls))
print(dim(pancreatitis))
print(dim(pancreatic_cancer))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (normal_controls)
# Assuming all matrices have the same number of features, use the number of rows from normal_controls
p <- dim(normal_controls)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 3  # Total number of classes
Y[[1]] <- t(normal_controls)
Y[[2]] <- t(pancreatitis)
Y[[3]] <- t(pancreatic_cancer)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 10 Artery Pressures
{# GSE24988_series_matrix.txt.gz
# Step 1: Read the entire file
file_content <- readLines("GSE24988_series_matrix.txt")

# Step 2: Extract metadata (Sample Titles and Characteristics)
sample_title_line <- grep("!Sample_title", file_content)
sample_titles <- unlist(strsplit(gsub("\"", "", file_content[sample_title_line]), "\t"))[-1]

sample_characteristics_line <- grep("!Sample_characteristics_ch2", file_content)
sample_characteristics <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_line]), "\t"))[-1]

# Step 3: Read the expression matrix
matrix_start <- grep("!series_matrix_table_begin", file_content)
matrix_data <- read.delim(textConnection(file_content[(matrix_start + 1):length(file_content)]), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 4: Add metadata as rows to the expression matrix
# Create a new dataframe with sample titles and characteristics as the first two rows
metadata_df <- data.frame(
  ID_REF = c("Sample_Title", "Sample_Characteristics"),
  rbind(sample_titles, sample_characteristics),
  stringsAsFactors = FALSE
)

# Ensure that the column names of the metadata match the matrix data (excluding ID_REF)
colnames(metadata_df) <- colnames(matrix_data)

# Step 5: Combine the metadata and the expression matrix
combined_data <- rbind(metadata_df, matrix_data)

dim(combined_data)

# Step 6: Save the combined data as a CSV file
# write.csv(combined_data, "GSE24988_combined_data.csv", row.names = FALSE)

# Step 6: Create matrices for each disease state group (excluding the last row)
# Remove the last row and convert values to numeric, handling NA and non-numeric data

severe_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("severe PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))
intermediate_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("intermediate PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))
no_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("no PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))
validation_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("validation set for PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))

# Step 7: Replace NA values with 0 in all matrices
severe_ph_matrix[is.na(severe_ph_matrix)] <- 0
intermediate_ph_matrix[is.na(intermediate_ph_matrix)] <- 0
no_ph_matrix[is.na(no_ph_matrix)] <- 0
validation_ph_matrix[is.na(validation_ph_matrix)] <- 0

# Step 8: Output the dimensions of each matrix to verify
print(dim(severe_ph_matrix))
print(dim(intermediate_ph_matrix))
print(dim(no_ph_matrix))
print(dim(validation_ph_matrix))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# GSE24988_series_matrix.txt.gz (4 classes)
p <- dim(severe_ph_matrix)[1]
Y <- list()
k <- 4
Y[[1]] <- t(severe_ph_matrix)
Y[[2]] <- t(intermediate_ph_matrix)
Y[[3]] <- t(no_ph_matrix)
Y[[4]] <- t(validation_ph_matrix)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3]),rep(0,sum(n)),rep(1,n[4])), ncol=k, nrow=sum(n))
C <- cbind(diag(q),-rep(1,q))


# normal-reference approach based test for GLHT problem
ZGZ2017.GLHT.2cNRT(Y,G,n,p)
ZZZ2022.GLHT.2cNRT(Y,X,C,n,p)
ZZG2022.GLHTBF.2cNRT(Y,G,n,p)
ZZ2022.GLHTBF.3cNRT(Y,G,n,p)
ZZ2022.GLHT.3cNRT(Y,G,n,p)

# Others' tests for GLHT problem
FHW2004.GLHT.NABT(Y,X,C,n,p)
SF2006.GLHT.NABT(Y,X,C,n,p)
S2007.ks.NABT(Y,n,p)
YS2012.GLHT.NABT(Y,X,C,n,p)
ZGZ2017.GLHTBF.NABT(Y,G,n,p)


# data 11 skeletal muscle
{# GSE25941_series_matrix.txt
# Load the Artery Pressures dataset (GSE24988)
# The dataset is in a compressed .gz file and is a series matrix format
# Step 1: Read the entire file
file_content <- readLines("GSE24988_series_matrix.txt")

# Step 2: Extract metadata (Sample Titles and Characteristics)
# Extract sample titles and characteristics to use as metadata for the dataset
sample_title_line <- grep("!Sample_title", file_content)
sample_titles <- unlist(strsplit(gsub("\"", "", file_content[sample_title_line]), "\t"))[-1]

sample_characteristics_line <- grep("!Sample_characteristics_ch2", file_content)
sample_characteristics <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_line]), "\t"))[-1]

# Step 3: Read the expression matrix
# Locate the start of the data matrix and load the matrix data from the file
matrix_start <- grep("!series_matrix_table_begin", file_content)
matrix_data <- read.delim(textConnection(file_content[(matrix_start + 1):length(file_content)]), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 4: Add metadata as rows to the expression matrix
# Create a new dataframe with sample titles and characteristics as the first two rows
metadata_df <- data.frame(
  ID_REF = c("Sample_Title", "Sample_Characteristics"),
  rbind(sample_titles, sample_characteristics),
  stringsAsFactors = FALSE
)

# Ensure that the column names of the metadata match the matrix data (excluding ID_REF)
colnames(metadata_df) <- colnames(matrix_data)

# Step 5: Combine the metadata and the expression matrix
combined_data <- rbind(metadata_df, matrix_data)

# Output the dimensions of the combined data to verify
dim(combined_data)

# Step 6: Create matrices for each disease state group (excluding the last row)
# Remove the last row and convert values to numeric, handling NA and non-numeric data
severe_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("severe PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))
intermediate_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("intermediate PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))
no_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("no PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))
validation_ph_matrix <- apply(combined_data[-c(1, 2, nrow(combined_data)), grep("validation set for PH", combined_data[2, ])], 2, function(x) as.numeric(as.character(x)))

# Step 7: Replace NA values with 0 in all matrices
severe_ph_matrix[is.na(severe_ph_matrix)] <- 0
intermediate_ph_matrix[is.na(intermediate_ph_matrix)] <- 0
no_ph_matrix[is.na(no_ph_matrix)] <- 0
validation_ph_matrix[is.na(validation_ph_matrix)] <- 0

# Output the dimensions of each matrix to verify correct data extraction
print(dim(severe_ph_matrix))
print(dim(intermediate_ph_matrix))
print(dim(no_ph_matrix))
print(dim(validation_ph_matrix))

#summary(young_females_matrix)
#summary(young_males_matrix)
#summary(old_females_matrix)
#summary(old_males_matrix)
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (severe_ph_matrix)
# Assuming all matrices have the same number of features, use the number of rows from severe_ph_matrix
p <- dim(severe_ph_matrix)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 4  # Total number of classes
Y[[1]] <- t(severe_ph_matrix)
Y[[2]] <- t(intermediate_ph_matrix)
Y[[3]] <- t(no_ph_matrix)
Y[[4]] <- t(validation_ph_matrix)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 12 SARS-CoV-2
{# GSE156063_series_matrix.txt
# Load the SARS-CoV-2 dataset (GSE156063)
# The dataset is split between a metadata file (.txt) and a gene expression data file (.csv)
# Step 1: Read the series matrix file (metadata file)
file_content <- readLines("GSE156063_series_matrix.txt")

# Step 2: Extract Sample Characteristics (Sample Characteristics ch2: disease states)
# Extract sample characteristics ch2 to identify disease states of the samples
sample_characteristics_ch2_line <- grep("!Sample_characteristics_ch2", file_content)

# Split the characteristics data into individual elements by tab separator
sample_characteristics_ch2 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch2_line]), "\t"))[-1]

# Step 3: Read the gene expression data from the CSV file (GSE156063_swab_gene_counts.csv)
# Load the gene expression data from the CSV file
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
# Uncomment the line below to save the combined data as a CSV file
# write.csv(combined_GSE156063, "GSE156063_combined_data.csv", row.names = FALSE)
dim(combined_GSE156063)

# Step 8: Create matrices for each disease state group (SC2, other viral, no virus)
# Select samples corresponding to each group based on disease state
sc2_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: SC2"])  # SC2 samples
other_virus_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: other virus"])  # Other viral samples
no_virus_matrix <- as.matrix(combined_GSE156063[-1, combined_GSE156063[1, ] == "disease state: no virus"])  # Non-viral samples

# Convert character data to numeric if necessary
sc2_matrix <- apply(sc2_matrix, 2, as.numeric)
other_virus_matrix <- apply(other_virus_matrix, 2, as.numeric)
no_virus_matrix <- apply(no_virus_matrix, 2, as.numeric)

# Step 9: Output the dimensions of each matrix to verify
print(dim(sc2_matrix))
print(dim(other_virus_matrix))
print(dim(no_virus_matrix))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (sc2_matrix)
# Assuming all matrices have the same number of features, use the number of rows from sc2_matrix
p <- dim(sc2_matrix)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 3  # Total number of classes
Y[[1]] <- t(sc2_matrix)
Y[[2]] <- t(other_virus_matrix)
Y[[3]] <- t(no_virus_matrix)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)



# data 13
{# GSE157103_series_matrix.txt
# Load the SARS-CoV-2 dataset (GSE157103)
# The dataset is split between a metadata file (.txt) and a gene expression data file (.tsv)
# Step 1: Read the metadata from the .txt file
file_content <- readLines("GSE157103_series_matrix.txt")

# Step 2: Extract metadata (Sample Characteristics Ch11 and Ch12)
# Extract sample characteristics ch11 and ch12 to identify the ICU status and disease state of the samples
sample_characteristics_ch11_line <- grep("!Sample_characteristics_ch11", file_content)
sample_characteristics_ch11 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch11_line]), "\t"))[-1]

sample_characteristics_ch12_line <- grep("!Sample_characteristics_ch12", file_content)
sample_characteristics_ch12 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch12_line]), "\t"))[-1]

# Step 3: Read the expression matrix from the .tsv file
# Load the gene expression data from the .tsv file
expression_matrix <- read.delim("GSE157103_genes.tpm.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 4: Add metadata as rows to the expression matrix
# Convert sample characteristics ch11 and ch12 to rows and combine with expression matrix
metadata_df <- data.frame(
  ID_REF = c("Sample_Characteristics_Ch11", "Sample_Characteristics_Ch12"),
  rbind(sample_characteristics_ch11, sample_characteristics_ch12),
  stringsAsFactors = FALSE
)

# Ensure that the column names of the metadata match the matrix data (excluding ID_REF)
colnames(metadata_df) <- colnames(expression_matrix)

# Step 5: Combine the metadata and the expression matrix
combined_GSE157103 <- rbind(metadata_df, expression_matrix)

# Step 6: Save the combined data as a CSV file (optional)
# Uncomment the line below to save the combined data as a CSV file
# write.csv(combined_GSE157103, "GSE157103_combined_data.csv", row.names = FALSE)
dim(combined_GSE157103)

# Step 7: Create matrices for each group based on the disease state and ICU status
# Ensure numeric coercion and handle potential non-numeric values

# Remove metadata rows (ensure only numeric data is used)
covid_nonicu_matrix <- as.matrix(
  apply(combined_GSE157103[-c(1, 2),
                           combined_GSE157103[1, ] == "disease state: COVID-19" &
                             combined_GSE157103[2, ] == "icu: no"], 2, as.numeric)
)

covid_icu_matrix <- as.matrix(
  apply(combined_GSE157103[-c(1, 2),
                           combined_GSE157103[1, ] == "disease state: COVID-19" &
                             combined_GSE157103[2, ] == "icu: yes"], 2, as.numeric)
)

noncovid_nonicu_matrix <- as.matrix(
  apply(combined_GSE157103[-c(1, 2),
                           combined_GSE157103[1, ] == "disease state: non-COVID-19" &
                             combined_GSE157103[2, ] == "icu: no"], 2, as.numeric)
)

noncovid_icu_matrix <- as.matrix(
  apply(combined_GSE157103[-c(1, 2),
                           combined_GSE157103[1, ] == "disease state: non-COVID-19" &
                             combined_GSE157103[2, ] == "icu: yes"], 2, as.numeric)
)

# Step 8: Output the dimensions of each matrix to verify
print(dim(covid_nonicu_matrix))
print(dim(covid_icu_matrix))
print(dim(noncovid_nonicu_matrix))
print(dim(noncovid_icu_matrix))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (covid_nonicu_matrix)
# Assuming all matrices have the same number of features, use the number of rows from covid_nonicu_matrix
p <- dim(covid_nonicu_matrix)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 4  # Total number of classes
Y[[1]] <- t(covid_nonicu_matrix)
Y[[2]] <- t(covid_icu_matrix)
Y[[3]] <- t(noncovid_nonicu_matrix)
Y[[4]] <- t(noncovid_icu_matrix)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# data 14
{ # rats.xls
# Load the rats dataset
# The dataset is in an Excel file (.xls) containing information on different rat groups
# Step 1: Load the necessary library to read the Excel file
library(readxl)

# Step 2: Read the excel file containing the rats data
# Load the data from the specified Excel file
rats_data <- read_excel("rats.xls")

# Step 3: Remove the GROUP and T_0 columns from the dataset
# Remove the columns that are not needed for hypothesis testing
rats_data_clean <- rats_data[, -c(1, 2)]

# Step 4: Split the data into groups based on the GROUP column values (1, 2, 3, 4)
# Create separate matrices for each group, transposing the data to match the expected format for hypothesis testing
group1_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 1, ]))
group2_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 2, ]))
group3_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 3, ]))
group4_matrix <- t(as.matrix(rats_data_clean[rats_data$GROUP == 4, ]))

# Step 5: Output the dimensions of each group matrix to verify that the split is correct
# Print the dimensions of each group matrix to ensure correct data extraction
print(dim(group1_matrix))
print(dim(group2_matrix))
print(dim(group3_matrix))
print(dim(group4_matrix))
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (group1_matrix)
# Assuming all matrices have the same number of features, use the number of rows from group1_matrix
p <- dim(group1_matrix)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 4  # Total number of classes
Y[[1]] <- t(group1_matrix)
Y[[2]] <- t(group2_matrix)
Y[[3]] <- t(group3_matrix)
Y[[4]] <- t(group4_matrix)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)


# Load the HDNRA library for hypothesis testing
# data 15
{ # eating
# Load the eating dataset from OpenML
# The dataset is retrieved from OpenML and split into different food classes for analysis
library(OpenML)
library(farff)

# Step 1: Load the dataset from OpenML
# Load the dataset using the specified data ID from OpenML
dataset <- getOMLDataSet(data.id = 1233)

# Step 2: Extract the data from the dataset object
# Extract the data frame from the dataset object
data <- dataset$data

# Step 3: Split the data based on the class label
# Ensure that 'class' is the correct column name for the labels
# Remove the class column from the data to focus on feature variables
class_column <- "class"  # Adjust if necessary based on the actual label column name
data_numeric <- data[, -which(names(data) == class_column)]  # Remove the class column from the data

# Split the data into different classes, ensuring only numeric data is used
class_list <- split(data_numeric, data[[class_column]])

# Step 4: Convert each class to a matrix and transpose
# Create transposed matrices for each class for hypothesis testing
matrix_apple <- t(as.matrix(class_list$Apple))
matrix_banana <- t(as.matrix(class_list$Banana))
matrix_biscuit <- t(as.matrix(class_list$Biscuit))
matrix_crisp <- t(as.matrix(class_list$Crisp))
matrix_haribo <- t(as.matrix(class_list$Haribo))
matrix_nectarine <- t(as.matrix(class_list$Nectarine))
matrix_no_food <- t(as.matrix(class_list$`No_Food`))  # Note the use of backticks for special characters

# Step 5: Check for NA values in each matrix
# Output the number of NA values in each matrix to verify data quality
sum(is.na(matrix_apple))
sum(is.na(matrix_banana))
sum(is.na(matrix_biscuit))
sum(is.na(matrix_crisp))
sum(is.na(matrix_haribo))
sum(is.na(matrix_nectarine))
sum(is.na(matrix_no_food))

# Step 6: Check the dimensions of each matrix
# Print the dimensions of each matrix to ensure correct data extraction
dim(matrix_apple)
dim(matrix_banana)
dim(matrix_biscuit)
dim(matrix_crisp)
dim(matrix_haribo)
dim(matrix_nectarine)
dim(matrix_no_food)
}

# Load the HDNRA library for hypothesis testing
library(HDNRA)

# Define the number of features (genes) based on the first matrix (matrix_apple)
# Assuming all matrices have the same number of features, use the number of rows from matrix_apple
p <- dim(matrix_apple)[1]

# Create a list of matrices, each containing the transposed feature data for each class
# Transposing to match the expected format for hypothesis testing (features as columns)
Y <- list()
k <- 7  # Total number of classes
Y[[1]] <- t(matrix_apple)
Y[[2]] <- t(matrix_banana)
Y[[3]] <- t(matrix_biscuit)
Y[[4]] <- t(matrix_crisp)
Y[[5]] <- t(matrix_haribo)
Y[[6]] <- t(matrix_nectarine)
Y[[7]] <- t(matrix_no_food)

# Define the number of samples for each class (number of rows in each transposed matrix)
n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]), nrow(Y[[5]]), nrow(Y[[6]]), nrow(Y[[7]]))

# Create the contrast matrix for group comparisons in General Linear Hypothesis Testing
# The matrix G allows for testing differences between groups, with the last group as reference
G <- cbind(diag(k-1), rep(-1, k-1))

# Set up parameters for General Linear Hypothesis Testing (GLHT)
q <- k - 1  # Define the number of contrasts

# Construct the design matrix for hypothesis testing (group membership for all samples)
# The matrix X defines group assignment for each sample
X <- matrix(c(rep(1, n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), 
              rep(1, n[3]), rep(0, sum(n)), rep(1, n[4]), rep(0, sum(n)), 
              rep(1, n[5]), rep(0, sum(n)), rep(1, n[6]), rep(0, sum(n)), 
              rep(1, n[7])), ncol=k, nrow=sum(n))

# Construct the contrast matrix for hypothesis testing
# C matrix represents linear combinations of model coefficients for testing specific hypotheses
C <- cbind(diag(q), -rep(1, q))

# Perform hypothesis tests using normal-reference approach-based methods
# Normal-reference approach based tests for GLHT problem
ZGZ2017.GLHT.2cNRT(Y, G, n, p)
ZZZ2022.GLHT.2cNRT(Y, X, C, n, p)
ZZG2022.GLHTBF.2cNRT(Y, G, n, p)
ZZ2022.GLHTBF.3cNRT(Y, G, n, p)
ZZ2022.GLHT.3cNRT(Y, G, n, p)

# Perform hypothesis tests using other methods for GLHT problem
# These functions provide alternative approaches for testing the General Linear Hypothesis
FHW2004.GLHT.NABT(Y, X, C, n, p)
SF2006.GLHT.NABT(Y, X, C, n, p)
S2007.ks.NABT(Y, n, p)
YS2012.GLHT.NABT(Y, X, C, n, p)
ZGZ2017.GLHTBF.NABT(Y, G, n, p)

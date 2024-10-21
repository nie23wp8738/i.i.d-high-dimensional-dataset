# empirical sizes for GLHT functions

# Load the Pancreatic dataset (GSE24279)
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

test_group <- t(pancreatic_cancer)
dim(test_group)

# Helper function to extract p.value and parameters from NRtest object
extract_results <- function(result) {
  # Check if the result is of class NRtest
  if (inherits(result, "NRtest")) {
    p.value <- result$p.value
    parameters <- result$parameter
    return(list(p.value = p.value, parameters = parameters))
  } else {
    return(list(p.value = NA, parameters = rep(NA, 2)))  # Return NA if the result is not valid
  }
}

# Load required libraries
library(doParallel)
library(HDNRA)
library(foreach)

# Set the significance level
alpha <- 0.10

# Set up parallel computing environment
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Number of repetitions for the analysis
nrep <- 10000

# Perform size control using 10 different GLHT functions with random groupings in each loop
results <- foreach(N = 1:nrep, .combine = 'rbind', .packages = c("HDNRA")) %dopar% {
  # Randomly split data into three groups for each iteration
  total_samples <- nrow(test_group)
  idx <- sample(1:total_samples)

  # Define split points for the groups
  split_1 <- floor(total_samples / 3)
  split_2 <- 2 * split_1
  p <- dim(test_group)[2]
  # Create three random groups
  Y <- list()
  Y[[1]] <- test_group[idx[1:split_1], ]
  Y[[2]] <- test_group[idx[(split_1+1):(split_2)], ]
  Y[[3]] <- test_group[idx[(split_2+1):total_samples], ]

  # Number of samples in each group
  n <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]))

  # Set up the design matrices for each method
  G <- cbind(diag(2), rep(-1, 2))  # Adjust for k = 3 groups
  q <- 2
  X <- matrix(c(rep(1,n[1]), rep(0, sum(n)), rep(1, n[2]), rep(0, sum(n)), rep(1, n[3])), ncol=3, nrow=sum(n))
  C <- cbind(diag(q), -rep(1, q))

  # Apply the 10 GLHT methods without error handling
  ZGZ2017_GLHT_2cNRT <- extract_results(ZGZ2017.GLHT.2cNRT(Y, G, n, p))
  ZZZ2022_GLHT_2cNRT <- extract_results(ZZZ2022.GLHT.2cNRT(Y, X, C, n, p))
  ZZG2022_GLHTBF_2cNRT <- extract_results(ZZG2022.GLHTBF.2cNRT(Y, G, n, p))
  ZZ2022_GLHTBF_3cNRT <- extract_results(ZZ2022.GLHTBF.3cNRT(Y, G, n, p))
  ZZ2022_GLHT_3cNRT <- extract_results(ZZ2022.GLHT.3cNRT(Y, G, n, p))
  FHW2004_GLHT_NABT <- extract_results(FHW2004.GLHT.NABT(Y, X, C, n, p))
  SF2006_GLHT_NABT <- extract_results(SF2006.GLHT.NABT(Y, X, C, n, p))
  S2007_ks_NABT <- extract_results(S2007.ks.NABT(Y, n, p))
  YS2012_GLHT_NABT <- extract_results(YS2012.GLHT.NABT(Y, X, C, n, p))
  ZGZ2017_GLHTBF_NABT <- extract_results(ZGZ2017.GLHTBF.NABT(Y, G, n, p))

  # Convert p-values to binary (0 or 1) based on significance level
  size_ZGZ2017_GLHT_2cNRT <- as.integer(ZGZ2017_GLHT_2cNRT$p.value < alpha)
  size_ZZZ2022_GLHT_2cNRT <- as.integer(ZZZ2022_GLHT_2cNRT$p.value < alpha)
  size_ZZG2022_GLHTBF_2cNRT <- as.integer(ZZG2022_GLHTBF_2cNRT$p.value < alpha)
  size_ZZ2022_GLHTBF_3cNRT <- as.integer(ZZ2022_GLHTBF_3cNRT$p.value < alpha)
  size_ZZ2022_GLHT_3cNRT <- as.integer(ZZ2022_GLHT_3cNRT$p.value < alpha)
  size_FHW2004_GLHT_NABT <- as.integer(FHW2004_GLHT_NABT$p.value < alpha)
  size_SF2006_GLHT_NABT <- as.integer(SF2006_GLHT_NABT$p.value < alpha)
  size_S2007_ks_NABT <- as.integer(S2007_ks_NABT$p.value < alpha)
  size_YS2012_GLHT_NABT <- as.integer(YS2012_GLHT_NABT$p.value < alpha)
  size_ZGZ2017_GLHTBF_NABT <- as.integer(ZGZ2017_GLHTBF_NABT$p.value < alpha)

  # Return sizes and parameters for all 10 tests
  return(c(
    size_ZGZ2017_GLHT_2cNRT, size_ZZZ2022_GLHT_2cNRT, size_ZZG2022_GLHTBF_2cNRT,
    size_ZZ2022_GLHTBF_3cNRT, size_ZZ2022_GLHT_3cNRT, size_FHW2004_GLHT_NABT,
    size_SF2006_GLHT_NABT, size_S2007_ks_NABT, size_YS2012_GLHT_NABT,
    size_ZGZ2017_GLHTBF_NABT,
    ZGZ2017_GLHT_2cNRT$parameters[1],
    ZZZ2022_GLHT_2cNRT$parameters,
    ZZG2022_GLHTBF_2cNRT$parameters[1],
    ZZ2022_GLHTBF_3cNRT$parameters[1],
    ZZ2022_GLHT_3cNRT$parameters[1],
    YS2012_GLHT_NABT$parameters
  ))
}

# Stop the cluster
stopCluster(cl)

# Convert results to matrices
results <- as.matrix(results)
size <- as.matrix(results[, 1:10])  # First 10 columns are the size results
para <- as.matrix(results[, 11:16]) # Next 10 columns are the parameters

# Print the proportion of significant results (size control)
print(colSums(size) / nrep)

# Print the average of parameters
print(colSums(para) / nrep)

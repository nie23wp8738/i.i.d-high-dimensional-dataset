# high-dimensional datasets

# data 1
{# yeoh-2002-v2_database.txt
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
}

library(HDNRA)

p <- dim(T_ALL)[1]
Y <- list()
k <- 6
Y[[1]] <- t(BCR_ABL)
Y[[2]] <- t(E2A_PBX1)
Y[[3]] <- t(Hyperdiploid_50)
Y[[4]] <- t(MLL)
Y[[5]] <- t(TEL_AML1)
Y[[6]] <- t(T_ALL)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]),nrow(Y[[6]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3]),rep(0,sum(n)),rep(1,n[4]),rep(0,sum(n)),rep(1,n[5]),rep(0,sum(n)),rep(1,n[6])),ncol=k,nrow=sum(n))
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


# data 2
{
# alizadeh-2000-v3_database.txt
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

# Check if there are any remaining NA values
# cat("DLBCL1 NA count:", sum(is.na(DLBCL1)), "\n")
# cat("DLBCL2 NA count:", sum(is.na(DLBCL2)), "\n")
# cat("FL NA count:", sum(is.na(FL)), "\n")
# cat("CLL NA count:", sum(is.na(CLL)), "\n")
}

library(HDNRA)

# alizadeh-2000-v3_database.txt (4 classes)
p <- dim(DLBCL1)[1]
Y <- list()
k <- 4
Y[[1]] <- t(DLBCL1)
Y[[2]] <- t(DLBCL2)
Y[[3]] <- t(FL)
Y[[4]] <- t(CLL)

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


# data 3
{
# tomlins-2006_database.txt
# Load the dataset and skip the first row
data3 <- read.table("tomlins-2006_database.txt", header=TRUE, sep="\t", fill=TRUE, skip=1)

# Verify the data structure
# head(data3)
# colnames(data3)

# Create matrices for each sample type
EPI <- as.matrix(data3[, grep("EPI", names(data3))])
MET <- as.matrix(data3[, grep("MET", names(data3))])
PCA <- as.matrix(data3[, grep("PCA", names(data3))])
PIN <- as.matrix(data3[, grep("PIN", names(data3))])
STROMA <- as.matrix(data3[, grep("STROMA", names(data3))])

# Output the dimensions of each matrix to verify
print(dim(EPI))
print(dim(MET))
print(dim(PCA))
print(dim(PIN))
print(dim(STROMA))
}

library(HDNRA)

# tomlins-2006_database.txt (5 classes)
p <- dim(EPI)[1]
Y <- list()
k <- 5
Y[[1]] <- t(EPI)
Y[[2]] <- t(MET)
Y[[3]] <- t(PCA)
Y[[4]] <- t(PIN)
Y[[5]] <- t(STROMA)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3]),rep(0,sum(n)),rep(1,n[4]),rep(0,sum(n)),rep(1,n[5])), ncol=k, nrow=sum(n))
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


# data 4
{
# Iris Data.csv
# Load the dataset
data_iris <- read.csv("Iris Data.csv")

# Transpose the dataset, setting species as the column names
data_iris_t <- t(data_iris[, 1:4])  # Transpose the features
colnames(data_iris_t) <- data_iris$species  # Set species as the column names

# Verify the transposed dataset
# head(data_iris_t)

# Create matrices for each species
Setosa <- as.matrix(data_iris_t[, grep("setosa", colnames(data_iris_t), ignore.case = TRUE)])
Versicolor <- as.matrix(data_iris_t[, grep("versicolor", colnames(data_iris_t), ignore.case = TRUE)])
Virginica <- as.matrix(data_iris_t[, grep("virginica", colnames(data_iris_t), ignore.case = TRUE)])

# Output the dimensions of each matrix to verify
print(dim(Setosa))
print(dim(Versicolor))
print(dim(Virginica))
}

library(HDNRA)

# Iris Data (3 classes)
p <- dim(Setosa)[1]
Y <- list()
k <- 3
Y[[1]] <- t(Setosa)
Y[[2]] <- t(Versicolor)
Y[[3]] <- t(Virginica)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3])), ncol=k, nrow=sum(n))
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


# data 5
{
# ecoli_data.csv
# Load the dataset
data_ecoli <- read.csv("ecoli_data.csv")

# Transpose the dataset, setting class as the column names
data_ecoli_t <- t(data_ecoli[, 1:7])  # Transpose the features
colnames(data_ecoli_t) <- data_ecoli$class  # Set class as the column names

# Verify the transposed dataset
# head(data_ecoli_t)

# Create matrices for each class
cp <- as.matrix(data_ecoli_t[, grep("cp", colnames(data_ecoli_t), ignore.case = TRUE)])
im <- as.matrix(data_ecoli_t[, grep("im", colnames(data_ecoli_t), ignore.case = TRUE)])  # Match "im" but not "imL", "imU", "imS"
pp <- as.matrix(data_ecoli_t[, grep("pp", colnames(data_ecoli_t), ignore.case = TRUE)])
imU <- as.matrix(data_ecoli_t[, grep("imU", colnames(data_ecoli_t), ignore.case = TRUE)])
om <- as.matrix(data_ecoli_t[, grep("om", colnames(data_ecoli_t), ignore.case = TRUE)])  # Match "om" but not "omL"
omL <- as.matrix(data_ecoli_t[, grep("omL", colnames(data_ecoli_t), ignore.case = TRUE)])
imS <- as.matrix(data_ecoli_t[, grep("imS", colnames(data_ecoli_t), ignore.case = TRUE)])
imL <- as.matrix(data_ecoli_t[, grep("imL", colnames(data_ecoli_t), ignore.case = TRUE)])

# Output the dimensions of each matrix to verify
print(dim(cp))
print(dim(im))
print(dim(pp))
print(dim(imU))
print(dim(om))
print(dim(omL))
print(dim(imS))
print(dim(imL))
}

library(HDNRA)

# ecoli_data.csv (5 classes)
p <- dim(cp)[1]
Y <- list()
k <- 5
Y[[1]] <- t(cp)
Y[[2]] <- t(im)
Y[[3]] <- t(pp)
Y[[4]] <- t(imU)
Y[[5]] <- t(om)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3]),rep(0,sum(n)),rep(1,n[4]),rep(0,sum(n)),rep(1,n[5])), ncol=k, nrow=sum(n))
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


# data 6
{
# yeast_data.csv
# Load the dataset
data_yeast <- read.csv("yeast_data.csv")

# Transpose the dataset, setting class as the column names
data_yeast_t <- t(data_yeast[, 1:8])  # Transpose the features
colnames(data_yeast_t) <- data_yeast$class  # Set class as the column names

# Verify the transposed dataset
# head(data_yeast_t)

# Create matrices for each class
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

# Output the dimensions of each matrix to verify
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

library(HDNRA)

# yeast_data.csv (9 classes)
p <- dim(CYT)[1]
Y <- list()
k <- 9
Y[[1]] <- t(CYT)
Y[[2]] <- t(NUC)
Y[[3]] <- t(MIT)
Y[[4]] <- t(ME3)
Y[[5]] <- t(ME2)
Y[[6]] <- t(ME1)
Y[[7]] <- t(EXC)
Y[[8]] <- t(VAC)
Y[[9]] <- t(POX)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]),nrow(Y[[6]]),nrow(Y[[7]]),nrow(Y[[8]]),nrow(Y[[9]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3]),rep(0,sum(n)),rep(1,n[4]),rep(0,sum(n)),rep(1,n[5]),rep(0,sum(n)),rep(1,n[6]),rep(0,sum(n)),rep(1,n[7]),rep(0,sum(n)),rep(1,n[8]),rep(0,sum(n)),rep(1,n[9])), ncol=k, nrow=sum(n))
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


# data 7 Heart Disease
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
}

library(HDNRA)

# processed_cleveland_data.csv (5 classes)
p <- dim(class_0)[1]
Y <- list()
k <- 5
Y[[1]] <- t(class_0)
Y[[2]] <- t(class_1)
Y[[3]] <- t(class_2)
Y[[4]] <- t(class_3)
Y[[5]] <- t(class_4)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3]),rep(0,sum(n)),rep(1,n[4]),rep(0,sum(n)),rep(1,n[5])), ncol=k, nrow=sum(n))
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


# data 8
{# wine.csv
# Load the dataset
data_wine <- read.csv("wine.csv")

# Transpose the dataset, setting the class as the column names
data_wine_t <- t(data_wine[, 2:ncol(data_wine)])  # Transpose the features, excluding the class column
colnames(data_wine_t) <- data_wine[, 1]  # Set the class column as the column names

# Verify the transposed dataset
# head(data_wine_t)

# Create matrices for each class
class_11 <- as.matrix(data_wine_t[, grep("1", colnames(data_wine_t))])
class_22 <- as.matrix(data_wine_t[, grep("2", colnames(data_wine_t))])
class_33 <- as.matrix(data_wine_t[, grep("3", colnames(data_wine_t))])

# Output the dimensions of each matrix to verify
print(dim(class_11))
print(dim(class_22))
print(dim(class_33))
}


library(HDNRA)

# wine.csv (3 classes)
p <- dim(class_11)[1]
Y <- list()
k <- 3
Y[[1]] <- t(class_11)
Y[[2]] <- t(class_22)
Y[[3]] <- t(class_33)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3])), ncol=k, nrow=sum(n))
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


# data 9 Pancreatic
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
}

library(HDNRA)

# GSE24279_non-normalized.txt.gz (3 classes)
p <- dim(normal_controls)[1]
Y <- list()
k <- 3
Y[[1]] <- t(normal_controls)
Y[[2]] <- t(pancreatitis)
Y[[3]] <- t(pancreatic_cancer)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3])), ncol=k, nrow=sum(n))
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

# Step 6: Save the combined data as a CSV file (optional)
# write.csv(combined_GSE25941, "GSE25941_combined_data.csv", row.names = FALSE)

# Step 7: Verify the result
dim(combined_GSE25941)
# head(combined_GSE25941)

# Step 8: Create matrices for each group (young females, young males, old females, old males)
# Exclude the last row using -nrow(combined_GSE25941)
young_females_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Female" & combined_GSE25941[1, ] == "age: Young"])
young_males_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Male" & combined_GSE25941[1, ] == "age: Young"])
old_females_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Female" & combined_GSE25941[1, ] == "age: Old"])
old_males_matrix <- as.matrix(combined_GSE25941[-c(1, 2, nrow(combined_GSE25941)), combined_GSE25941[2, ] == "gender: Male" & combined_GSE25941[1, ] == "age: Old"])

young_females_matrix <- apply(young_females_matrix, 2, function(x) as.numeric(as.character(x)))
young_males_matrix <- apply(young_males_matrix, 2, function(x) as.numeric(as.character(x)))
old_females_matrix <- apply(old_females_matrix, 2, function(x) as.numeric(as.character(x)))
old_males_matrix <- apply(old_males_matrix, 2, function(x) as.numeric(as.character(x)))

#sum(is.na(young_females_matrix))
#sum(is.na(young_males_matrix))
#sum(is.na(old_females_matrix))
#sum(is.na(old_males_matrix))


# Step 9: Standardize (scale) the matrices for each group
# This ensures that each feature has a mean of 0 and a standard deviation of 1
young_females_matrix <- scale(young_females_matrix)
young_males_matrix <- scale(young_males_matrix)
old_females_matrix <- scale(old_females_matrix)
old_males_matrix <- scale(old_males_matrix)

# Step 10: Output the dimensions of each matrix to verify correctness
print(dim(young_females_matrix))
print(dim(young_males_matrix))
print(dim(old_females_matrix))
print(dim(old_males_matrix))

#summary(young_females_matrix)
#summary(young_males_matrix)
#summary(old_females_matrix)
#summary(old_males_matrix)
}

library(HDNRA)

# GSE25941_series_matrix.txt (4 classes)
p <- dim(young_females_matrix)[1]
Y <- list()
k <- 4
Y[[1]] <- t(young_females_matrix)
Y[[2]] <- t(young_males_matrix)
Y[[3]] <- t(old_females_matrix)
Y[[4]] <- t(old_males_matrix)

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


# data 12 SARS-CoV-2
{# GSE156063_series_matrix.txt
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


library(HDNRA)

# GSE156063_series_matrix.txt (3 classes)
p <- dim(sc2_matrix)[1]
Y <- list()
k <- 3
Y[[1]] <- t(sc2_matrix)
Y[[2]] <- t(other_virus_matrix)
Y[[3]] <- t(no_virus_matrix)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3])), ncol=k, nrow=sum(n))
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


# data 13
{# GSE157103_series_matrix.txt
# GSE157103_series_matrix.txt and GSE157103_genes.tpm.tsv
# Step 1: Read the metadata from the .txt file
file_content <- readLines("GSE157103_series_matrix.txt")

# Step 2: Extract metadata (Sample Characteristics Ch11 and Ch12)
sample_characteristics_ch11_line <- grep("!Sample_characteristics_ch11", file_content)
sample_characteristics_ch11 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch11_line]), "\t"))[-1]

sample_characteristics_ch12_line <- grep("!Sample_characteristics_ch12", file_content)
sample_characteristics_ch12 <- unlist(strsplit(gsub("\"", "", file_content[sample_characteristics_ch12_line]), "\t"))[-1]

# Step 3: Read the expression matrix from the .tsv file
expression_matrix <- read.delim("GSE157103_genes.tpm.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#expression_matrix <- expression_matrix[, -1]  # Exclude the first column containing gene symbols

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

library(HDNRA)

# GSE157103_series_matrix.txt (4 classes)
p <- dim(covid_nonicu_matrix)[1]
Y <- list()
k <- 4
Y[[1]] <- t(covid_nonicu_matrix)
Y[[2]] <- t(covid_icu_matrix)
Y[[3]] <- t(noncovid_nonicu_matrix)
Y[[4]] <- t(noncovid_icu_matrix)

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


# data 14
{ # rats.xls
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
}


library(HDNRA)

# rats.xls (4 classes)
p <- dim(group1_matrix)[1]
Y <- list()
k <- 4
Y[[1]] <- t(group1_matrix)
Y[[2]] <- t(group2_matrix)
Y[[3]] <- t(group3_matrix)
Y[[4]] <- t(group4_matrix)

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



# data 15
{ # eating
library(OpenML)
library(farff)
# Step 1: Load the dataset from OpenML
dataset <- getOMLDataSet(data.id = 1233)

# Step 2: Extract the data from the dataset object
data <- dataset$data

# Step 3: Split the data based on the class label
# Ensure that 'class' is the correct column name for the labels
class_column <- "class"  # Adjust if necessary based on the actual label column name
data_numeric <- data[, -which(names(data) == class_column)]  # Remove the class column from the data

# Split the data into different classes, ensuring only numeric data is used
class_list <- split(data_numeric, data[[class_column]])

# Step 4: Convert each class to a matrix and transpose
matrix_apple <- t(as.matrix(class_list$Apple))
matrix_banana <- t(as.matrix(class_list$Banana))
matrix_biscuit <- t(as.matrix(class_list$Biscuit))
matrix_crisp <- t(as.matrix(class_list$Crisp))
matrix_haribo <- t(as.matrix(class_list$Haribo))
matrix_nectarine <- t(as.matrix(class_list$Nectarine))
matrix_no_food <- t(as.matrix(class_list$`No_Food`))  # Note the use of backticks for special characters


# Step 6: Check for NA values in each matrix
sum(is.na(matrix_apple))
sum(is.na(matrix_banana))
sum(is.na(matrix_biscuit))
sum(is.na(matrix_crisp))
sum(is.na(matrix_haribo))
sum(is.na(matrix_nectarine))
sum(is.na(matrix_no_food))

# Step 7: Check the dimensions of each matrix
dim(matrix_apple)
dim(matrix_banana)
dim(matrix_biscuit)
dim(matrix_crisp)
dim(matrix_haribo)
dim(matrix_nectarine)
dim(matrix_no_food)
}

library(HDNRA)

# eating dataset (7 classes)
p <- dim(matrix_apple)[1]
Y <- list()
k <- 7
Y[[1]] <- t(matrix_apple)
Y[[2]] <- t(matrix_banana)
Y[[3]] <- t(matrix_biscuit)
Y[[4]] <- t(matrix_crisp)
Y[[5]] <- t(matrix_haribo)
Y[[6]] <- t(matrix_nectarine)
Y[[7]] <- t(matrix_no_food)

n <- c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]),nrow(Y[[6]]),nrow(Y[[7]]))
G <- cbind(diag(k-1),rep(-1,k-1))

q <- k-1
X <- matrix(c(rep(1,n[1]),rep(0,sum(n)),rep(1,n[2]),rep(0,sum(n)),rep(1,n[3]),rep(0,sum(n)),rep(1,n[4]),rep(0,sum(n)),rep(1,n[5]),rep(0,sum(n)),rep(1,n[6]),rep(0,sum(n)),rep(1,n[7])), ncol=k, nrow=sum(n))
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

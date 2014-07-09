# Big data chain ladder with activeH5

# Load the package
require(rhdf5)
require(ChainLadder)


#--------------------------------------------------------------------------------------------------
# Creating the chain ladder data
#--------------------------------------------------------------------------------------------------

# Age-To-Age Factors
ageFact <- seq(1.9, 1, by = -.1)

# Inflation Rate
infRate <- 1.02

revCols <- function(x){
  x[,ncol(x):1]
}

# We use shake rather than jitter, shake is much faster and equivalent to my python Jitter() function
shake <- function(vec, sigmaScale = 100)
{
  rnorm(n = length(vec), mean = vec, sd = vec/sigmaScale)
}

# Alternative Row generation funtion
GenerateRow <- function(iDev, dFactors = cumprod(ageFact), dInflationRate = 1.02, initClaim = 154)
{
  shake(initClaim)*shake(c(1, dFactors))*(dInflationRate^iDev)
}

# Function to generate a claims matrix
GenerateTriangle <- function(iSize, ...)
{
  indices = 1:iSize
  mClaimTri = t(sapply(indices, GenerateRow, ...))
  # Reverse columns to get the claims triangle
  mClaimTri = revCols(mClaimTri)
  # Assign nan to lower triangle
  mClaimTri[lower.tri(mClaimTri)] = NA
  mClaimTri = revCols(mClaimTri)
  return(mClaimTri)
}


# Creating the claims H5 file
claimsH5 <- "data/ClaimsTri.h5"
claimsH5File <- H5Fcreate(claimsH5)

# Function to write a triangle to H5 file
WriteToH5File <- function(sName, h5File = claimsH5File){
  h5writeDataset.matrix(GenerateTriangle(11), h5File, name = sName, level = 0)
}

sMatrixNames <- as.character(1:2000)
system.time(null <- lapply(sMatrixNames, WriteToH5File))
#  user  system elapsed 
# 7.773   0.025   7.795 

#--------------------------------------------------------------------------------------------------
# Analysing triangles using the chain ladder package
#--------------------------------------------------------------------------------------------------

# H5 File to store the claims triangle
sProcessedH5 <- "data/ClaimsSquare.h5"
processedH5 <- H5Fcreate(sProcessedH5)

# Function to process an item in the ChainLadder file
ChainLadder <- function(sName = "1", file = claimsH5File){
  MackChainLadder(h5read(file = file, name = sName), est.sigma="Mack")$FullTriangle
}

# Function to write a processed chainladder to H5 file
WriteSquare <- function(sName = "1", file = processedH5){
  h5writeDataset.matrix(ChainLadder(sName = sName), processedH5, name = sName, level = 0)
}

system.time(null <- lapply(sMatrixNames, WriteSquare))
#  user  system elapsed 
# 203.940 813.666 131.732

#--------------------------------------------------------------------------------------------------
# Analysing the chain ladder data
#--------------------------------------------------------------------------------------------------

# Function for calculating age-to-age factors
GetFactor <- function(index, mTri)
{
  fact = matrix(mTri[-c((nrow(mTri) - index + 1):nrow(mTri)), index:(index + 1)], ncol = 2)
  fact = c(sum(fact[,1]), sum(fact[,2]))  
  return(fact[2]/fact[1])
}

GetChainSquare <- function(mClaimTri)
{
  nCols <- ncol(mClaimTri)
  dFactors = sapply(1:(nCols - 1), GetFactor, mTri = mClaimTri)
  dAntiDiag = diag(revCols(mClaimTri))[2:nCols]
  for(index in 1:length(dAntiDiag))
    mClaimTri[index + 1, (nCols - index + 1):nCols] = dAntiDiag[index]*cumprod(dFactors[(nCols - index):(nCols - 1)])
  mClaimTri
}

# Alternative technique
ChainLadder2 <- function(sName = "1", file = claimsH5File){
  GetChainSquare(h5read(file = file, name = sName))
}


# Function to write a processed chainladder to H5 file
WriteSquare2 <- function(sName = "1", file = processedH5){
  h5writeDataset.matrix(ChainLadder2(sName = sName), processedH5, name = sName, level = 0)
}

# H5 File to store the claims triangle
sProcessedH5 <- "data/ClaimsSquare.h5"
processedH5 <- H5Fcreate(sProcessedH5)
system.time(null <- lapply(sMatrixNames, WriteSquare2))
#  user  system elapsed 
#  12.957   0.047  12.986 


# Closing the files
H5Fclose(claimsH5File)
H5Fclose(processedH5)

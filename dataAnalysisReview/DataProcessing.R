# An R file for processing data.
# We consider three different datasets.

# 1. Breast_GSE70947
# 2. Liver_GSE14520
# 3. Prostate_GSE6919

# For data processing we have the following goals:
#   1. Load the dataset
#   2. Find the indices for abnormal/cancer samples
#   3. Recode factor to binary 0/1, where 1 is abnormal sample
#   4. Return data with first column as response others as predictors.
process.dat <- function(name = "Breast_GSE70947_small") {
  # Drop the first column.
  dat <- data.table::fread(file = paste0("data/",name, ".csv"), drop = 1)
  n <- nrow(dat)
  names(dat)[1] <- "y"

  ind_norm <- which(dat$y =="normal")

  dat$y[ind_norm] <- 0
  dat$y[-ind_norm] <- 1
  dat$y <- as.factor(dat$y)

  # Remove this gene for breast cancer data due to
  # computational issues with SpAM
  if(name == "Breast_GSE70947") {
    dat[,"NM_001195605"] <- NULL
  }
  return(dat)
}


# A data processing for other datasets
process.dat2 <- function(name = "spam") {
  # Drop the first column.
  if(name == "spam") {
    data(spam, package = "kernlab")
  }
  dat <- spam[,-58]
  dat <- cbind(spam$type, dat)
  n <- nrow(dat)
  names(dat)[1] <- "y"

  ind_norm <- which(dat$y =="nonspam")

  dat$y <- as.numeric(dat$y)
  dat$y[ind_norm] <- 0
  dat$y[-ind_norm] <- 1
  dat$y <- as.factor(dat$y)

  return(dat)
}

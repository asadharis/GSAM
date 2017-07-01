# This file loads the data for analysis.

# Load clinical data
clinical <- read.csv("data/clinical.csv", na.strings = "[Not Available]")
rna <- read.csv("data/RNA.csv", check.names = FALSE)
rna.gene <- rna$Gene
rna.id <- rna$EntrezID
rna <- rna[,-(1:2)]

# Generate the ids of subjects from the RNA data
string.list <- strsplit(colnames(rna), split = "-")
string.list <- string.list[-(1:2)]

rna.ids <- sapply(string.list, function(x){
  paste0(x[1],"-",x[2],"-",x[3])
  })

# We now find the intersection of the two list of names
inter.ids <- intersect(rna.ids, clinical$bcr_patient_barcode)

ind.clinical <- which(clinical$bcr_patient_barcode %in% inter.ids)
ind.rna <- which(rna.ids %in% inter.ids)

# Generate the response vector
y <- clinical$psa_result_preop[ind.clinical]
y <- clinical$psa_value[ind.clinical]# $psa_result_preop[ind.clinical]
# Obtain covariate matrix
x <- rna[,ind.rna]
x <- as.data.frame(t(x))
names(x) <- rna.gene

dat <- as.data.frame(cbind("PSA" = y, x))

# Finally we obtain the data set with complete cases.
dat <- dat[complete.cases(dat),]

#inds <- which(apply(dat,2,function(x){all(x>0)}))
#dat <- dat[,inds]

dat <- log(dat + 0.1)

save("dat", file = "DataFile.RData")


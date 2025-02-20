# script analyzing manipulate gene expression data
# setwd("D:/Rstudio gene expression analysis/R script")
setwd("D:/Rstudio gene expression analysis/R script")

#installing packages RNA-seq Data analysis
install.packages("dplyr")
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
# Install cluster package to user-specific library
install.packages("cluster", lib = "C:/Users/LENOVO/AppData/Local/R/win-library/4.4")

# Verify installation
if ("GEOquery" %in% rownames(installed.packages())) {
  print("GEOquery is installed and loaded successfully.")
} else {
  print("GEOquery is not installed.")
}

# load libraries
library("dplyr")
library("tidyverse")
library("GEOquery")

# Install conflicted package if not already installed
install.packages("conflicted")

# Load conflicted package
library("conflicted")

# Set conflict preferences
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")

# read in the data
dat <- read.csv(file ="D:/Rstudio gene expression analysis/R script/GSE183947_fpkm.csv")
dim(dat)

#get metadata
gse <- getGEO (GEO = 'GSE183947', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# Check the column names, then deleted it 
colnames(metadata)
metadata$`tissue:ch1`<- NULL

#select specific column
metadata.subset <- select(metadata, c(1, 10, 11, 17))
nrow(metadata.subset)
head(metadata.subset)

#rename column rename function and edit specific variables mutate function
metadata.modified <- metadata %>%
  select (1, 10, 11, 17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename (metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue:", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis:", "", metastasis)) %>%
  head(60)

head(metadata.modified)

head(dat)

#reshaping data
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)
head(dat.long)

#join dataframes = dat.long + metadata.modified
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))
head(dat.long)

#explore data
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)

  

  



  
  








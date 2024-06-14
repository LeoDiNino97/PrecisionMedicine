# 0 - Load packages --------------

BiocManager::install("BiocGenerics")
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
install.packages("GGally")

required_packages <- c("BiocGenerics", "DESeq2", "psych", "NetworkToolbox", "ggplot2",
                       "GGally", "sna", "network", "TCGAbiolinks", 
                       "SummarizedExperiment", "DT", "latex2exp", "gridExtra")

lapply(required_packages, library, character.only = TRUE)

seed <- 123
set.seed(seed)

colors <- c("#2EBFA5", "#F26419", "#FDCA40")


# 1 - Data loading and preprocessing -----------------------------------------------------------

# Project ID
proj <- "TCGA-LIHC"
directory <- proj

# Look for all data linked to a "Primary Tumor" sample and store them in a dataframe
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")

GDCdownload(query = rna.query.C, directory = directory, method = "api")
rna.data.C <- GDCprepare(rna.query.C, directory = directory)
rna.expr.data.C <- assay(rna.data.C)
genes.info.c <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))


# Apply the same procedure to the "Normal tissue"
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

GDCdownload(query = rna.query.N, directory = directory, method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = directory)
rna.expr.data.N <- assay(rna.data.N)
genes.info.n <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))

# Retrieving patients data
clinical.query <- GDCquery_clinic(project = proj, 
                                  type = "clinical", 
                                  save.csv = FALSE)

# Before moving on, let's look if genes.info.c and genes.info.n are both needed
all.equal(genes.info.c, genes.info.n)

# Since they are duplicates, we keep just one of them to release memory

genes.info <- genes.info.c
rm(genes.info.n)
rm(genes.info.c)


# 2 - Data normalization with DeSeq2 ------------------------------------------

all(rownames(rna.expr.data.C) == rownames(rna.expr.data.N))
full.data <- cbind(rna.expr.data.C, rna.expr.data.N)

dim(full.data)
dim(rna.expr.data.C)
dim(rna.expr.data.N)

full.data <- data.frame(full.data)

# Defining meta data for the normalization procedure
metad <- rep("normal", dim(full.data)[2])
metad[1:dim(rna.expr.data.C)[2]] <- "cancer"
metad <- data.frame(metad)

rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad

metad[,1] <- as.factor(metad[,1])
full.data <- cbind(rownames(full.data), full.data)

dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=metad, 
                              design= ~ condition,
                              tidy=TRUE)

View(counts(dds))
dim(counts(dds))

# Perform normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, 
                            normalized=TRUE)

sum(rowSums(normalized_counts == 0) == dim(full.data)[2]) # No null rows :)

# Split back the table
rna.expr.data.C <- as.data.frame(normalized_counts[, 1:dim(rna.expr.data.C)[2]])
rna.expr.data.N <- as.data.frame(normalized_counts[, (dim(rna.expr.data.C)[2]+1):dim(normalized_counts)[2]])


# 3 - Data engineering and cleaning -----------------------------------------------------------

# Let's read the patients which we are interested into
patients <- read.table("annex_1.txt", stringsAsFactors=FALSE)[,1]

# Looking for duplicates
dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) # No duplicates :)
dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) # No duplicates :)


# Renaming patients to match the given format for patient_id
colnames(rna.expr.data.C) <- gsub("\\.", "-", substr(colnames(rna.expr.data.C), 1,12))
unique(colnames(rna.expr.data.C))

colnames(rna.expr.data.N) <- gsub("\\.", "-", substr(colnames(rna.expr.data.N), 1,12))
unique(colnames(rna.expr.data.N))

# Let's match the two datasets with the given patients_id list
rna.expr.data.C <- rna.expr.data.C[,match(patients, colnames(rna.expr.data.C))] 
rna.expr.data.N <- rna.expr.data.N[,match(patients, colnames(rna.expr.data.N))]

# Patients id of patients being in both tables
length(intersect(colnames(rna.expr.data.C), colnames(rna.expr.data.N))) == 49

# All given patients have both cancer and normal tissues data

# Let's retrieve gender based identifiers
women <- intersect(patients,
                   clinical.query[clinical.query$gender == "female",]$submitter_id)

men <- intersect(patients,
                 clinical.query[clinical.query$gender == "male",]$submitter_id)

expr.C.women <- rna.expr.data.C[,women]
expr.C.men <- rna.expr.data.C[,men]

expr.N.women <- rna.expr.data.N[,women]
expr.N.men <- rna.expr.data.N[,men]

# Select the genes of interest

genes <- read.table("gene-list2.txt", stringsAsFactors=FALSE)[,1]

genes.id <- genes.info[genes.info$gene_name %in% genes, "gene_id"]

genes.c <- intersect(rownames(rna.expr.data.C), 
                     genes.id) 

genes.n <- intersect(rownames(rna.expr.data.N),  
                     genes.id)

# Filter out our the previously retrieved tables 

expr.C.men <- expr.C.men[genes.c,]
expr.C.women <- expr.C.women[genes.c,]

expr.N.men <- expr.N.men[genes.n,]
expr.N.women <- expr.N.women[genes.n,]


# Differentially Expressed Genes (DEGs) -----------------------------------



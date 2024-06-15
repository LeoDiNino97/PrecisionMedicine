# 0 - Load packages --------------
#install.packages("psych")
#install.packages("BiocManager")
#install.packages("NetworkToolbox")
#install.packages("sna")
library(BiocManager)
library(psych)
library(NetworkToolbox)
library(sna)
library(DT)

BiocManager::install("BiocGenerics")
BiocManager::install("TCGAbiolinks")

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
dir.create(file.path(proj))

# Look for all data linked to a "Primary Tumor" sample and store them in a dataframe
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")

GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")
rna.expr.data.C <- assay(rna.data.C)
genes.info.c <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))


# Apply the same procedure to the "Normal tissue"
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata" )
rna.expr.data.N <- assay(rna.data.N)
genes.info.n <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))

all(na.omit(genes.info.n) == na.omit(genes.info.c))

clinical.query<- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

View(clinical.query)
table(clinical.query$ajcc_pathologic_stage)

boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query,
        col = "gold", main = "Title", xlab = "", ylab= "age", las=2 )

# Before moving on, let's look if genes.info.c and genes.info.n are both needed
all.equal(genes.info.c, genes.info.n)

# Since they are duplicates, we keep just one of them to release memory

genes.info <- genes.info.c
rm(genes.info.n)
rm(genes.info.c)


# Data normalization with DeSeq2 ------------------------------------------

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


# Data engineering and cleaning -----------------------------------------------------------
setwd("C:\\Users\\GiuseppeDiPoce\\Desktop\\DEPM")
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

# Get rid of ENSG00000197976.12_PAR_Y 
row <- "ENSG00000197976.12_PAR_Y"

expr.N.men <- expr.N.men[!rownames(expr.N.men) %in% row,]
expr.N.women <- expr.N.men[!rownames(expr.N.women) %in% row,]
expr.C.men <- expr.N.men[!rownames(expr.C.men) %in% row,]
expr.C.women <- expr.N.men[!rownames(expr.C.women) %in% row,]
dim(expr.C.women)

# Differentially Expressed Genes (DEGs) -----------------------------------

# Analysis for men data

# Fold change computations
fc.men <- log2(rowMeans(expr.C.men) / rowMeans(expr.N.men) ) 
names(fc.men) <- rownames(expr.C.men)

row <- "ENSG00000197976.12_PAR_Y"
fc.men <- fc.men[!names(fc.men) %in% row]

# t-statistics p-value computation and multiple comparisons adjustments
pval.fc.men <- sapply(1:nrow(expr.C.men), function(i) (t.test(expr.C.men[i,], 
                                                              expr.N.men[i,]))$p.value)
pval.fc.men.fdr <- p.adjust(pval.fc.men, method="fdr")
names(pval.fc.men.fdr) <- rownames(expr.C.men)
pval.fc.men.fdr <- pval.fc.men.fdr[!names(pval.fc.men.fdr) %in% row]

expr.table.men <- data.frame(cbind(fc.men, pval.fc.men.fdr))
expr.table.men[,1] <- round(expr.table.men[,1],2)

deg.genes.men <- rownames(expr.table.men[abs(expr.table.men$fc) >= 1.2 & expr.table.men$pval.fc.men.fdr <= 0.05,]) 
deg.genes.men

head(expr.table.men[deg.genes.men,], 10)

# Volcano plot

expr.table.men$diffexpressed <- "NO";
expr.table.men$diffexpressed[expr.table.men$fc >= 1.2 & expr.table.men$pval.fc.men.fdr <= 0.05] <- "UP"
expr.table.men$diffexpressed[expr.table.men$fc <= -1.2 & expr.table.men$pval.fc.men.fdr <= 0.05] <- "DOWN"
head(expr.table.men)

expr.table.men$diffexpressed <- as.factor(expr.table.men$diffexpressed)
summary(expr.table.men$diffexpressed)

p1 <- ggplot(data=expr.table.men, aes(x=fc.men, 
                                      y=-log10(pval.fc.men.fdr), 
                                      col=diffexpressed))+  
  geom_point(size=1.1) +
  theme_bw() +
  scale_color_manual(values = c("darkred", "grey", "darkblue")) +
  xlab("fold change (log2)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(0.05), col="black", lty='dotted')+
  geom_vline(xintercept=1.2,col="black", lty='dotted')+
  geom_vline(xintercept=-1.2, col="black", lty='dotted') + 
  labs(color = 'DEGs status ')

p1
cat(deg.genes.men, sep = "\n")

# Analysis for women data

# Fold change computations
fc.women <- log2(rowMeans(expr.C.women) / rowMeans(expr.N.women) ) 
names(fc.women) <- rownames(expr.C.women)
head(fc.women)

row <- "ENSG00000197976.12_PAR_Y"
fc.women <- fc.women[!names(fc.women) %in% row]

# t-statistics p-value computation and multiple comparision adjustments
pval.fc.women <- sapply(1:nrow(expr.C.women), function(i) (t.test(expr.C.women[i,], 
                                                                  expr.N.women[i,]))$p.value)
pval.fc.women.fdr <- p.adjust(pval.fc.women, method="fdr")
names(pval.fc.women.fdr) <- rownames(expr.C.women)
pval.fc.women.fdr <- pval.fc.women.fdr[!names(pval.fc.women.fdr) %in% row]

expr.table.women <- data.frame(cbind(fc.women, pval.fc.women.fdr))
expr.table.women[,1] <- round(expr.table.women[,1],2)

deg.genes.women <- rownames(expr.table.women[abs(expr.table.women$fc) >= 1.2 & expr.table.women$pval.fc.women.fdr <= 0.05,]) 
deg.genes.women

head(expr.table.men[deg.genes.women,], 10)

# Volcano plot

expr.table.women$diffexpressed <- "NO";
expr.table.women$diffexpressed[expr.table.women$fc >= 1.2 & expr.table.women$pval.fc.women.fdr <= 0.05] <- "UP"
expr.table.women$diffexpressed[expr.table.women$fc <= -1.2 & expr.table.women$pval.fc.women.fdr <= 0.05] <- "DOWN"
head(expr.table.women)

expr.table.women$diffexpressed <- as.factor(expr.table.women$diffexpressed)
summary(expr.table.women$diffexpressed)

p2 <- ggplot(data=expr.table.women, aes(x=fc.women, 
                                        y=-log10(pval.fc.women.fdr), 
                                        col=diffexpressed))+  
  geom_point(size=1.1) +
  theme_bw() +
  scale_color_manual(values = c("darkred", "grey", "darkblue")) +
  xlab("fold change (log2)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(0.05), col="black", lty='dotted')+
  geom_vline(xintercept=1.2,col="black", lty='dotted')+
  geom_vline(xintercept=-1.2, col="black", lty='dotted') + 
  labs(color = 'DEGs status ')

p2
cat(deg.genes.women, sep = "\n")

# For export 

legend <- get_legend(p1)

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

grid.arrange(
  arrangeGrob(p1, p2, ncol = 2),
  legend,
  ncol = 2,
  widths = c(4, 1)
)

overlapping.DEGs <- intersect(deg.genes.men, deg.genes.women)



# Patient Similarity Network ----------------------------------------------

# Rebuild back a full dataframe with men and women 

expr.C <- cbind(expr.C.men, expr.C.women)

# Correlation for the cancer network using Pearson correlation
cor.mat.cp <- corr.test(expr.C, 
                        use = "pairwise", 
                        method = "pearson", 
                        adjust="fdr", 
                        ci=FALSE)

# Matrix with the pearson correlation coefficient
rho.cp <- cor.mat.cp$r

# Set the diagonal to zero (we are not interested to the edge with a gene with itself)
diag(rho.cp) <- 0

# Binary matrix, keep the edge for big correlation
qval.cp <- cor.mat.cp$p
qval.cp[lower.tri(qval.cp)] <- t(qval.cp)[lower.tri(qval.cp)]

# Correlation network of cancer samples

# Since the pvalues are all extremely small we use the correlation matrix directly to build the graph
adj.mat.cp <- rho.cp * (abs(rho.cp) > 0.7)   


net.cp <- as.network(adj.mat.cp, directed = FALSE)

# Let's ensure that our graph is connected
l.comp.p <- component.largest(net.cp, result = "graph")
l.comp.p <- adj.mat.cp[rownames(l.comp.p), rownames(l.comp.p)]

write.csv2(l.comp.p, "input-matrix-c.csv")
# python btc-community-c.py input-matrix-c.csv

# Read results 

comm.res.p <- read.csv2("output.txt", header = FALSE)
rownames(comm.res.p) <- rownames(l.comp.p)
sort(table(comm.res.p[,1]), decreasing = T)
length(table(comm.res.p[,1]))

net.final.p <- network(l.comp.p, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)


# Assign communities
net.final.p  %v% "community" <-  as.character(comm.res.p[,1])

ncom <- length(unique(net.final.p  %v% "community"))


set.seed(420)
pal <- sample(colors(distinct = T), ncom)
names(pal) <- 1:ncom

node_mapping <- data.frame(real_label = unique(net.final.p %v% "vertex.names"),
                           new_label = 1:length(unique(net.final.p %v% "vertex.names")))

# Update the network with new labels (to better identify patients we label them with numbers from 1 to 18)
net.final.p %v% "vertex.names" <- node_mapping$new_label


# Plot subnetwork

ggnet2(net.final.p, color = "community", palette =  pal, alpha = 1, 
       size = 5, edge.color = "grey", edge.alpha = 1, edge.size = 0.15, label = TRUE, label.size = 5)+
  guides(size = "none") 






# Differential Co-expressed Network ---------------------------------------
#3.1 Computation. Compute the differential co-expression network

# Now we would like to build the co-expression network, to do so we first need
# to compute the adjacency matrix

# We will define the matrix using the correlation between the genes so
# the edge between two genes will be present if the Z (t-value associated with the correlation coefficient) 
# t is used to test the hypothesis that the correlation is significantly different from zero

Z <- 3.5# we choose to set a higher value of Z statistic due to a very dense behavior of the graph connectivity with lower values

# Before computing correlations, we log transform the data using log2(x+1)
#We will treat the two conditions separately men/woman
filtr.expr.c.men <- log2(expr.C.men+1)
filtr.expr.n.men <- log2(expr.N.men+1) 
filtr.expr.c.women <- log2(expr.C.women+1)
filtr.expr.n.women <- log2(expr.N.women+1) 

# Create the [CANCER NETWORK]
# Correlation matrix for the cancer network
cor.mat.cancer_men <- corr.test(t(filtr.expr.c.men), use = "pairwise", 
                       method = "pearson", adjust="fdr", ci=FALSE)
cor.mat.cancer_women <- corr.test(t(filtr.expr.c.women), use = "pairwise", 
                                method = "pearson", adjust="fdr", ci=FALSE)

# matrix with the Pearson correlation coefficient
Z.c_men <- cor.mat.cancer_men$t
Z.c_women <- cor.mat.cancer_women$t
# Set the diagonal to zero (no self loops)
diag(Z.c_men) <- 0
diag(Z.c_women) <- 0
# Binary matrix, keep the edge wrt Z stat
adj.mat.c_men <- 1 * (abs(Z.c_men) > Z)
adj.mat.c_women <- 1 * (abs(Z.c_women) > Z)

if (!require("igraph")) {
  install.packages("igraph", dependencies = TRUE)
  library(igraph)
}



# Same process for the NORMAL NETWORK
# Correlation for the normal network
cor.mat.n.men <- corr.test(t(filtr.expr.n.men), use = "pairwise", 
                       method = "pearson", adjust="fdr", ci=FALSE)
cor.mat.n.women <- corr.test(t(filtr.expr.n.women), use = "pairwise", 
                       method = "pearson", adjust="fdr", ci=FALSE)

Z.n.women <- cor.mat.n.women$t
Z.n.men <- cor.mat.n.men$t
diag(Z.n.men) <- 0
diag(Z.n.women) <- 0

adj.mat.n_women <- 1 * (abs(Z.n.women) > Z)
adj.mat.n_men <- 1 * (abs(Z.n.men) > Z)


#-------------------------------------------------------------------------------
#Create network for each genre and condition;

# Create the NORMAL NETWORK given the adjacency matrix
net.n_men <- as.network(adj.mat.n_men, directed = FALSE)
net.n_women <- as.network(adj.mat.n_women, directed=FALSE)
# Let's check the number of edges and nodes
network.size(net.n_men) #648
network.edgecount(net.n_men) # 14195
# and the density of the network
network.density(net.n_men)   #0.06771519


# Same for the CANCER NETWORK
# Create the NORMAL NETWORK given the adjacency matrix
net.c_men <- as.network(adj.mat.c_men, directed = FALSE)
net.c_women <- as.network(adj.mat.c_women, directed=FALSE)
# Let's check the number of edges and nodes
network.size(net.c_men) #648
network.edgecount(net.c_men)  #14191
# and the density of the network
network.density(net.c_men) #0.06769611

# The two share apparently same density of the network 

#Scale free : A network is considered "scale-free" when its degree distribution
#follows a power law, at least asymptotically. This means that in such networks,
#a few nodes (often called hubs) have a significantly higher number of 
#connections compared to most other nodes, which have very few links.

# First for the NORMAL graph: 
# extract the degree of each node
d.n <- sna::degree(net.n_men, gmode = 'graph')
names(d.n) <- network.vertex.names(net.n_men)

# Print the histogram of the degree together with a line for the 95% quantile
x_n <- quantile(d.n[d.n>0],0.95)
hist(d.n,col = "lightblue", main = "Degree distribution - normal tissue")
abline(v=x_n, col="red")
# It doesn't look scale free
# The number of hubs of the network is
hubs.n <- d.n[d.n>=x_n]
length(hubs.n) # 31

# Power law: Another way to check if a network is scale free is to verify if it
# follows the power law. i.e. there should be a link between the fraction of 
# nodes f(k) with degree k and k itself f(k) ⁓ k^(-γ) (usually 2< γ <3)
d.n.table <- table(d.n)

# Convert the table to a data frame
d.n.fd <- data.frame(degree = as.numeric(names(d.n.table)),
                     count = as.numeric(d.n.table)/length(hubs.n))

plot(log(d.n.fd$degree), log(d.n.fd$count), 
     xlab = "Degree", ylab = "Degree Count", 
     main = "Degree Distribution", 
     pch = 16, col = "lightblue")
# For sure not scale free because there is not linear relationship among the two axes



# [CANCER NETWORK]: 
# extract the degree of each node
d.n_c <- sna::degree(net.c_men, gmode = 'graph')
names(d.n_c) <- network.vertex.names(net.c_men)

# Print the histogram of the degree together with a line for the 95% quantile
x_n <- quantile(d.n[d.n>0],0.95)
hist(d.n_c,col = "lightblue", main = "Degree distribution - normal tissue")
abline(v=x_n, col="red")
# It doesn't look scale free
# The number of hubs of the network is
hubs.n_c <- d.n_c[d.n_c>=x_n]
length(hubs.n_c) # 33

# Power law: Another way to check if a network is scale free is to verify if it
# follows the power law. i.e. there should be a link between the fraction of 
# nodes f(k) with degree k and k itself f(k) ⁓ k^(-γ) (usually 2< γ <3)
d.n_c.table <- table(d.n_c)

# Convert the table to a data frame
d.n_c.fd <- data.frame(degree = as.numeric(names(d.n_c.table)),
                     count = as.numeric(d.n_c.table)/length(hubs.n_c))

plot(log(d.n_c.fd$degree), log(d.n_c.fd$count), 
     xlab = "Degree", ylab = "Degree Count", 
     main = "Degree Distribution", 
     pch = 16, col = "lightblue")
# For sure not scale free because there is not linear relationship among the two axes




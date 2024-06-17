# Environment setting --------------

BiocManager::install("BiocGenerics")
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
install.packages("GGally")

required_packages <- c("BiocGenerics", "DESeq2", "psych", "NetworkToolbox", 
                       "ggplot2", "GGally", "sna", "network", "TCGAbiolinks", 
                       "igraph", "SummarizedExperiment", "DT", "latex2exp", "gridExtra", 
                       "cowplot", "poweRlaw","intergraph")

lapply(required_packages, library, character.only = TRUE)

seed <- 123
set.seed(seed)

# Data loading and preprocessing -----------------------------------------------------------

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


# Data normalization within populations + filtering with genes of interest -----

all(rownames(rna.expr.data.C) == rownames(rna.expr.data.N))
full.data <- cbind(rna.expr.data.C, rna.expr.data.N)

dim(full.data)
dim(rna.expr.data.C)
dim(rna.expr.data.N)

full.data <- data.frame(full.data)

# Let's read the patients which we are interested into
patients <- read.table("annex_1.txt", stringsAsFactors=FALSE)[,1]

# Looking for duplicates
dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) # No duplicates :)
dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) # No duplicates :)

# patients.sub <- gsub("-", "\\.", patients)
# Renaming patients to match the given format for patient_id
colnames(rna.expr.data.C) <- gsub("\\.", "-", substr(colnames(rna.expr.data.C), 1,12))
unique(colnames(rna.expr.data.C))

colnames(rna.expr.data.N) <- gsub("\\.", "-", substr(colnames(rna.expr.data.N), 1,12))
unique(colnames(rna.expr.data.N))

# Let's match the two datasets with the given patients_id list
rna.expr.data.C <- rna.expr.data.C[,match(patients, substr(colnames(rna.expr.data.C), 1, 12))] 
rna.expr.data.N <- rna.expr.data.N[,match(patients, substr(colnames(rna.expr.data.N), 1, 12))]

# Patients id of patients being in both tables
length(intersect(
  substr(colnames(rna.expr.data.C),1,12), 
  substr(colnames(rna.expr.data.N),1,12))) == 49

# All given patients have both cancer and normal tissues data

# Let's retrieve gender based identifiers
women <- intersect(patients,
                   clinical.query[clinical.query$gender == "female",]$submitter_id)

men <- intersect(patients,
                 clinical.query[clinical.query$gender == "male",]$submitter_id)

# Split data according to gender
expr.C.women <- rna.expr.data.C[,women]
expr.C.men <- rna.expr.data.C[,men]

expr.N.women <- rna.expr.data.N[,women]
expr.N.men <- rna.expr.data.N[,men]


# Defining meta data for the normalization procedure

# Normalizing men data
full.men <- cbind(expr.C.men, expr.N.men)
colnames(full.men) <- gsub("-", "\\.", colnames(full.men))
colnames(full.men)[28:54] <- paste(colnames(full.men)[28:54], ".1", sep = "")

metad.men <- rep("normal", dim(full.men)[2])
metad.men[1:dim(expr.C.men)[2]] <- "cancer"
metad.men <- data.frame(metad.men)

rownames(metad.men) <- colnames(full.men)
colnames(metad.men)[1] <- "condition"

metad.men[,1] <- as.factor(metad.men[,1])
#full.men <- cbind(rownames(full.men), full.men)

dds.men <- DESeqDataSetFromMatrix(countData=full.men, 
                                  colData=metad.men, 
                                  design= ~ condition,
                                  tidy=FALSE)

# Perform normalization
dds.men <- estimateSizeFactors(dds.men)
normalized.men <- counts(dds.men, 
                            normalized=TRUE)


# Split back the table
expr.C.men <- as.data.frame(normalized.men[, 1:dim(expr.C.men)[2]])
expr.N.men <- as.data.frame(normalized.men[, (dim(expr.C.men)[2]+1):dim(normalized.men)[2]])

# Normalizing women data
full.women <- cbind(expr.C.women, expr.N.women)
colnames(full.women) <- gsub("-", "\\.", colnames(full.women))
colnames(full.women)[23:44] <- paste(colnames(full.women)[23:44], ".1", sep = "")

metad.women <- rep("normal", dim(full.women)[2])
metad.women[1:dim(expr.C.women)[2]] <- "cancer"
metad.women <- data.frame(metad.women)

rownames(metad.women) <- colnames(full.women)
colnames(metad.women)[1] <- "condition"

metad.women[,1] <- as.factor(metad.women[,1])
#full.men <- cbind(rownames(full.men), full.men)

dds.women <- DESeqDataSetFromMatrix(countData=full.women, 
                                  colData=metad.women, 
                                  design= ~ condition,
                                  tidy=FALSE)

# Perform normalization
dds.women <- estimateSizeFactors(dds.women)
normalized.women <- counts(dds.women, 
                          normalized=TRUE)


# Split back the table
expr.C.women <- as.data.frame(normalized.women[, 1:dim(expr.C.women)[2]])
expr.N.women <- as.data.frame(normalized.women[, (dim(expr.C.women)[2]+1):dim(normalized.women)[2]])

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
overlapping.DEGs.info <- genes.info[overlapping.DEGs,]
write.csv2(overlapping.DEGs.info, "overlappingDEGS.csv")

sum(expr.table.men[overlapping.DEGs,]$diffexpressed == "UP")

# Differential co-expression network  -------------------------------------

Z <- 2

# We will treat the two conditions separately men/woman
expr.C.men <- log2(expr.C.men+1)
expr.C.women <- log2(expr.C.women+1)

expr.N.men <- log2(expr.N.men+1)
expr.N.women <- log2(expr.N.women+1)

# Create the nets
# Correlation matrix for the cancer network

# Cancer wrt sex
cor.mat.C.men <- corr.test(t(expr.C.men), use = "pairwise", 
                                method = "pearson", adjust="fdr", ci=FALSE)
cor.mat.C.women <- corr.test(t(expr.C.women), use = "pairwise", 
                                  method = "pearson", adjust="fdr", ci=FALSE)

# Normal wrt sex
cor.mat.N.men <- corr.test(t(expr.N.men), use = "pairwise", 
                           method = "pearson", adjust="fdr", ci=FALSE)
cor.mat.N.women <- corr.test(t(expr.N.women), use = "pairwise", 
                             method = "pearson", adjust="fdr", ci=FALSE)

# Fisher transform function

Fischer.Z <- function(M) {
  return(0.5 * log((1+M)/(1-M)))
}

# Apply Fisher transform function to correlation matrices
cor.mat.C.men <- Fischer.Z(cor.mat.C.men$r)
cor.mat.C.women <- Fischer.Z(cor.mat.C.women$r)
cor.mat.N.men <- Fischer.Z(cor.mat.N.men$r)
cor.mat.N.women <- Fischer.Z(cor.mat.N.women$r)

diag(cor.mat.C.men) <- 0
diag(cor.mat.C.women) <- 0
diag(cor.mat.N.men) <- 0
diag(cor.mat.N.women) <- 0

# Compute Z-scores

# Men graph
# Sample size of men between condition C and condition N is the same by design
n.men <- dim(expr.C.men)[2]
Z.men <- (cor.mat.C.men - cor.mat.N.men)/(2*sqrt(1/(n.men - 3)))

adj.mat.men <- 1 * (abs(Z.men) > Z)

# Women graph
# Sample size of men between condition C and condition N is the same by design
n.women <- dim(expr.C.women)[2]
Z.women <- (cor.mat.C.women - cor.mat.N.women)/(2*sqrt(1/(n.women - 3)))

adj.mat.women <- 1 * (abs(Z.women) > Z)

# Build the differential expression networks

DEnet.men <- as.network(component.largest(as.network(adj.mat.men, directed = FALSE), result = "graph"))
DEnet.women <- as.network(component.largest(as.network(adj.mat.women, directed=FALSE), result = "graph"))

# Men analysis
network.size(DEnet.men) 
network.edgecount(DEnet.men) 
network.density(DEnet.men)

d.men <- sna::degree(DEnet.men, gmode = 'graph')
names(d.men) <- network.vertex.names(DEnet.men)

# Print the histogram of the degree together with a line for the 95% quantile
q.men <- quantile(d.men[d.men>0],0.95)
hist(d.men,col = "lightblue", breaks = 80, main="")
abline(v=q.men, col="red", lty = 'dotted')

hubs.men <- d.men[d.men>=q.men]
length(hubs.men) 

# Deeper into degree distribution
d.men.table <- table(d.men)

# Convert the table to a data frame
d.men.fd <- data.frame(degree = as.numeric(names(d.men.table)),
                       count = as.numeric(d.men.table))

plot(log(d.men.fd$degree), log(d.men.fd$count), 
     xlab = "log.Degree", ylab = "log.DegreeCount", 
     pch = 16, col = "lightblue")

x.men <- log(d.men.fd$degree)
y.men <- log(d.men.fd$count)

model.lm.men <- glm(y.men[2:length(y.men)] ~ x.men[2:length(x.men)])
model.lm.men$coefficients

abline(a = as.numeric(model.lm.men$coefficients[1]), 
       b = as.numeric(model.lm.men$coefficients[2]), 
       col = 'red', 
       lty='dotted')

# Women analysis
network.size(DEnet.women) 
network.edgecount(DEnet.women) 
network.density(DEnet.women) 

d.women <- sna::degree(DEnet.women, gmode = 'graph')
names(d.women) <- network.vertex.names(DEnet.women)

# Print the histogram of the degree together with a line for the 95% quantile
q.women <- quantile(d.women[d.women>0],0.95)
hist(d.women,col = "pink", breaks = 80, main = "")
abline(v=q.women, col="red", lty = 'dotted')

hubs.women <- d.women[d.women>=q.women]
length(hubs.women) 

# Deeper into degree distribution
d.women.table <- table(d.women)

# Convert the table to a data frame
d.women.fd <- data.frame(degree = as.numeric(names(d.women.table)),
                         count = as.numeric(d.women.table))

plot(log(d.women.fd$degree), log(d.women.fd$count), 
     xlab = "log.Degree", ylab = "log.DegreeCount", 
     pch = 16, col = "pink")

x.women <- log(d.women.fd$degree)
y.women <- log(d.women.fd$count)

model.lm.women <- glm(y.women[2:length(y.women)] ~ x.women[2:length(x.women)])
model.lm.women$coefficients

abline(a = as.numeric(model.lm.women$coefficients[1]), 
       b = as.numeric(model.lm.women$coefficients[2]), 
       col = 'red', 
       lty='dotted')

# Overlapping hubs
overlapping.hubs <- intersect(names(hubs.men), names(hubs.women))
overlapping.hubs

# Comparing fitted line with MLE for power law

pl.men <- displ$new(rep(d.men.fd$degree, d.men.fd$count))
est.men <- estimate_xmin(pl.men)
pl.men$setXmin(est.men)

alpha.men <- pl.men$pars
xmin.men <- pl.men$xmin
  
pl.women <- displ$new(rep(d.women.fd$degree, d.women.fd$count))
est.women <- estimate_xmin(pl.women)
pl.women$setXmin(est.women)

# Get the estimated parameters
alpha.women <- pl.women$pars
xmin.women <- pl.women$xmin

#_______________________________________________________________________________

# Female vs male : only for cancer condition DEnetwork
#compute corr matrices
cor.mat.C.men <- corr.test(t(expr.C.men), use = "pairwise", 
                           method = "pearson", adjust="fdr", ci=FALSE)
cor.mat.C.women <- corr.test(t(expr.C.women), use = "pairwise", 
                             method = "pearson", adjust="fdr", ci=FALSE)

cor.mat.C.men <- Fischer.Z(cor.mat.C.men$r)
cor.mat.C.women <- Fischer.Z(cor.mat.C.women$r)

# Compute Z-scores

# Sample size of men and women onto Cancer tissue condition
n.C_men <- dim(expr.C.men)[2]
n.C_women <- dim(expr.C.women)[2]
Z.C <- (cor.mat.C.men - cor.mat.C.women)/(sqrt(1/(n.C_men - 3)) + (sqrt(1/(n.C_women - 3)) ) )

adj.mat.C <- 1 * (abs(Z.C) > Z)
table(adj.mat.C)


DEnet.C <- as.network(component.largest(as.network(adj.mat.C, directed=FALSE), result="graph"))

network.size(DEnet.C) 
network.edgecount(DEnet.C) 
network.density(DEnet.C) 

d.C <- sna::degree(DEnet.C, gmode = 'graph')
names(d.C) <- network.vertex.names(DEnet.C)

# Print the histogram of the degree together with a line for the 95% quantile
q.C <- quantile(d.C[d.C>0],0.95)
hist(d.C,col = "lightgreen", breaks=80, main = "")
abline(v=q.C, col="red", lty = 'dotted')

hubs.C <- d.C[d.C>=q.C]
length(hubs.C) # 9

#let's check differences among overlapping hubs compute before and those computed above 
intersect(names(hubs.C), names(overlapping.hubs)) #None :(


d.C.table <- table(d.C)

# Convert the table to a data frame
d.C.fd <- data.frame(degree = as.numeric(names(d.C.table)),
                       count = as.numeric(d.C.table))

plot(log(d.C.fd$degree), log(d.C.fd$count), 
     xlab = "log.Degree", ylab = "log.DegreeCount", 
     main = "", 
     pch = 16, col = "lightgreen")

x.C <- log(d.C.fd$degree)
y.C <- log(d.C.fd$count)

model.lm.C <- glm(y.C[2:length(y.C)] ~ x.C[2:length(x.C)])
model.lm.C$coefficients
abline(a = as.numeric(model.lm.C$coefficients[1]), 
       b = as.numeric(model.lm.C$coefficients[2]), col = 'red', lty='dotted')

intersect(names(hubs.men), hubs.C)
intersect(names(hubs.women), hubs.C)

# MLE for power-law
pl.C <- displ$new(rep(d.C.fd$degree, d.C.fd$count))
est.C <- estimate_xmin(pl.C)
pl.C$setXmin(est.C)

alpha.C <- pl$pars
xmin.C <- pl$xmin


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
adj.mat.cp <- rho.cp * (abs(rho.cp) > 0.85)   


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
net.final.p  %v% "Community" <-  as.character(comm.res.p[,1])

ncom <- length(unique(net.final.p  %v% "Community"))


set.seed(420)
pal <- sample(colors(distinct = T), ncom)
names(pal) <- 1:ncom

node_mapping <- data.frame(real_label = unique(net.final.p %v% "vertex.names"),
                           new_label = 1:length(unique(net.final.p %v% "vertex.names")))

net.final.p %v% "vertex.names" <- node_mapping$new_label


# Plot subnetwork

ggnet2(net.final.p, 
       color = "Community", 
       palette = pal, 
       alpha = 10, 
       size = 5, 
       edge.color = "lightgrey", 
       edge.alpha = 1, 
       edge.size = 0.15, 
       label = TRUE, 
       label.size = 3) 

# Some insights on the retrieved communities
node_mapping$community <- as.character(comm.res.p[,1])
node_mapping$real_label <- gsub("\\.","-",node_mapping$real_label)
communities <- clinical.query[clinical.query$submitter_id %in% node_mapping$real_label,]
community.1 <- communities[communities$submitter_id %in% node_mapping[node_mapping$community == 1,]$real_label,]
community.2 <- communities[communities$submitter_id %in% node_mapping[node_mapping$community == 2,]$real_label,]

prop.table(table(community.1$ajcc_pathologic_stage))
prop.table(table(community.2$ajcc_pathologic_stage))

prop.table(table(community.1$vital_status))
prop.table(table(community.2$vital_status))

prop.table(table(community.1$gender))
prop.table(table(community.2$gender))

prop.table(table(community.1$race))
prop.table(table(community.2$race))

prop.table(table(community.1$days_to_death))
prop.table(table(community.2$days_to_death))



# Optional tasks: differential expression networks for normal tissue ----------------------------------------------------------

cor.mat.N.men <- corr.test(t(expr.N.men), use = "pairwise", 
                           method = "pearson", adjust="fdr", ci=FALSE)
cor.mat.N.women <- corr.test(t(expr.N.women), use = "pairwise", 
                             method = "pearson", adjust="fdr", ci=FALSE)

cor.mat.N.men <- Fischer.Z(cor.mat.N.men$r)
cor.mat.N.women <- Fischer.Z(cor.mat.N.women$r)

# Compute Z-scores

# Sample size of men and women onto Cancer tissue condition
n.N.men <- dim(expr.N.men)[2]
n.N.women <- dim(expr.N.women)[2]
Z.N <- (cor.mat.N.men - cor.mat.N.women)/(sqrt(1/(n.N.men - 3)) + (sqrt(1/(n.N.women - 3)) ) )

adj.mat.N <- 1 * (abs(Z.N) > Z) #Try to change Z values and see how change the connectivity 
table(adj.mat.N)


DEnet.N <- as.network(component.largest(as.network(adj.mat.N, directed=FALSE), result='graph'))

# Network analysis
network.size(DEnet.N) 
network.edgecount(DEnet.N) 
network.density(DEnet.N)

d.N <- sna::degree(DEnet.N, gmode = 'graph')
names(d.N) <- network.vertex.names(DEnet.N)

# Print the histogram of the degree together with a line for the 95% quantile
q.N <- quantile(d.N[d.N>0],0.95)
hist(d.N,col = "gold", breaks=50, main = "")
abline(v=q.N, col="red", lty = 'dotted')

hubs.N <- d.N[d.N>=q.N]
length(hubs.N) 

d.N.table <- table(d.N)

# Convert the table to a data frame
d.N.fd <- data.frame(degree = as.numeric(names(d.N.table)),
                     count = as.numeric(d.N.table))


x.N <- log(d.N.fd$degree)
y.N <- log(d.N.fd$count)
x.N <- cbind(1, x.N[2:length(x.N)])
model.lm.N <- glm.fit(x.N, y.N[2:length(y.N)])
beta.n <- model.lm.N$coefficients
plot(x.N[,2], y.N[2:length(y.N)], 
     xlab = "log.Degree", ylab = "log.Count", main = "",pch = 16, col = "gold")

# Add the fitted line
abline(a = beta.n[1], b = beta.n[2], col = 'red', lty = 'dotted')

# MLE for power-law
pl.N <- displ$new(rep(d.N.fd$degree, d.N.fd$count))
est.N <- estimate_xmin(pl.N)
pl.N$setXmin(est.N)

alpha.N <- pl.N$pars
xmin.N <- pl.N$xmin

# Comparing with previous results

intersect(names(hubs.N), names(hubs.men))
intersect(names(hubs.N), names(hubs.women))

# Optional tasks: signed network analysis -------------------------------------------------

cor.mat.C.men <- corr.test(t(expr.C.men), use = "pairwise", 
                           method = "pearson", adjust="fdr", ci=FALSE)
cor.mat.C.women <- corr.test(t(expr.C.women), use = "pairwise", 
                             method = "pearson", adjust="fdr", ci=FALSE)

#fisher transform
cor.mat.C.men <- Fischer.Z(cor.mat.C.men$r)
cor.mat.C.women <- Fischer.Z(cor.mat.C.women$r)

n.C.men <- dim(expr.C.men)[2]
n.C.women <- dim(expr.C.women)[2]
Z.C <- (cor.mat.C.men - cor.mat.C.women)/(sqrt(1/(n.C.men - 3)) + (sqrt(1/(n.C.women - 3)) ) )
adj.C.pos <- 1 * (Z.C > Z) 
adj.C.neg <- 1 * (Z.C < -Z)
table(adj.C.neg)

#get the networks
DEnet.pos <- as.network(adj.C.pos, directed=FALSE)
DEnet.neg <- as.network(adj.C.neg, directed=FALSE)

# Plot netwroks
plot.network(DEnet.pos, displaylabels = FALSE, main = "Cancer Tissue positive Network")
plot.network(DEnet.neg, displaylabels = FALSE, main = "Cancer Tissue negative Network")

# Retrieve degree distribution
d.pos <- sna::degree(DEnet.pos, gmode = 'graph')
d.neg <- sna::degree(DEnet.neg, gmode='graph')
names(d.pos) <- network.vertex.names(DEnet.pos)
names(d.neg) <- network.vertex.names(DEnet.neg)

# Extract hubs for the two subnetworks
q_pos <- quantile(d.pos[d.pos>0],0.95)
hubs.pos <- names(d.pos)[d.pos>=q_pos]
length(hubs.pos) 


q_neg <- quantile(d.neg[d.neg>0],0.95)
hubs.neg <- names(d.neg)[d.neg>=q_neg]
length(hubs.neg) 

intersect(hubs.pos,hubs.neg)


# Subnetwork induced by hubs visualization --------------------------------

DEnet.men <- as.network(adj.mat.men, directed = FALSE)
DEnet.women <- as.network(adj.mat.women, directed=FALSE)
DEnet.C <- as.network(adj.mat.C, directed = FALSE)
DEnet.N <- as.network(adj.mat.N, directed=FALSE)

clean_matrix <- function(mat) {
  mat[is.na(mat) | is.nan(mat) | is.infinite(mat)] <- 0
  return(mat)
}

# Clean the adjacency matrices
adj.mat.men <- clean_matrix(adj.mat.men)
adj.mat.women <- clean_matrix(adj.mat.women)
adj.mat.C <- clean_matrix(adj.mat.C)
adj.mat.N <- clean_matrix(adj.mat.N)

g.men <- graph_from_adjacency_matrix(adj.mat.men, mode = "undirected")
g.women <- graph_from_adjacency_matrix(adj.mat.women, mode = "undirected")
g.N <- graph_from_adjacency_matrix(adj.mat.C, mode = "undirected")
g.C <- graph_from_adjacency_matrix(adj.mat.N, mode = "undirected")

neighbors.men <- unique(unlist(sapply(as.numeric(hubs.men), neighbors, graph = g.men)))
neighbors.women <- unique(unlist(sapply(as.numeric(hubs.women), neighbors, graph = g.women)))
neighbors.C <- unique(unlist(sapply(as.numeric(hubs.C), neighbors, graph = g.C)))
neighbors.N <- unique(unlist(sapply(as.numeric(hubs.N), neighbors, graph = g.N)))

neighbors.men <- unique(c(neighbors.men, as.numeric(hubs.men)))
neighbors.women <- unique(c(neighbors.women, as.numeric(hubs.women)))
neighbors.C <- unique(c(neighbors.C, as.numeric(hubs.C)))
neighbors.N <- unique(c(neighbors.N, as.numeric(hubs.N)))


subnet.plot <- function(g, suppnet, nodes, color) {
  # Create the edge data frame
  edge_df <- as_edgelist(g, names = TRUE)
  
  # Define colors and sizes
  node_colors <- ifelse(network.vertex.names(suppnet) %in% network.vertex.names(suppnet)[nodes], color, 'grey90')
  node_sizes <- ifelse(network.vertex.names(suppnet) %in% network.vertex.names(suppnet)[nodes], 1, 1)
  edge_colors <- ifelse(edge_df[,1] %in% network.vertex.names(suppnet)[nodes] & 
                        edge_df[,2] %in% network.vertex.names(suppnet)[nodes], color, "grey90")
  edge_sizes <- ifelse(edge_df[,1] %in% network.vertex.names(suppnet)[nodes] & 
                       edge_df[,2] %in% network.vertex.names(suppnet)[nodes], 1, 0.1)
  
  # Plot the network
  ggnet2(
    g,
    node.size = node_sizes,
    node.color = node_colors,
    edge.color = edge_colors,
    edge.size = edge_sizes,
  ) + 
    coord_equal() 
}


subnet.plot(g.men, DEnet.men, neighbors.men, 'darkblue')


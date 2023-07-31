# Loads required packages.
library(DESeq2)
library(apeglm)

# Use command to install packages if not installed.
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")

# Imports counts files and phenotype files for regular CLIP samples and combined CLIP/SM samples. 
# All counts files should have been edited as described on the doseCLIP Git page.
counts_clip <- as.matrix(read.csv("clip_subread_derived_counts.csv", row.names="gene_id"))
phenotype_clip <- read.csv("clip_phenotype_file.csv", row.names=1)
counts_clip_sm <- as.matrix(read.csv("clip_sm_subread_derived_counts.csv", row.names="gene_id"))
phenotype_clip_sm <- read.csv("clip_phenotype_file.csv", row.names=1)

# Imports counts file and imports file for RNA-Seq data. All counts files should have been 
# edited as described on the doseCLIP Git page.
counts_rna <- as.matrix(read.csv("rna_derived_counts.csv", row.names="gene_id"))
phenotype_rna <- read.csv("rna_phenotype_file.csv", row.names=1)

# Checks to make sure the sample names in the counts file matches the phenotype file.
all(rownames(counts_clip) == colnames(phenotype_clip))
all(rownames(counts_clip_sm) == colnames(phenotype_clip_sm))
all(rownames(counts_rna) == colnames(phenotype_rna))    

# Makes DeSeq2 object using counts and phenotype files.
dds_clip <- DESeqDataSetFromMatrix(countData = counts_clip, 
                                   colData = phenotype_clip, design = ~ condition)
dds_clip_sm <- DESeqDataSetFromMatrix(countData = counts_clip_sm, 
                                      colData = phenotype_clip_sm, design = ~ condition)
dds_rna <- DESeqDataSetFromMatrix(countData = counts_rna,
                                  colData = phenotype_rna, design = ~ condition)

# Filters out events with less than 20 reads per row and then performs differential
# gene expression analysis.
# CLIP.
keep <- rowSums(counts(dds_clip)) >= 20
dds_clip <- dds_clip[keep,]
dds_clip <- DESeq(dds_clip)
# CLIP/SM.
keep <- rowSums(counts(dds_clip_sm)) >= 20
dds_clip_sm <- dds_clip_sm[keep,]
dds_clip_sm <- DESeq(dds_clip_sm)
# RNA-Seq
keep <- rowSums(counts(dds_rna)) >= 20
dds_rna <- dds_rna[keep,]
dds_clip_sm <- DESeq(dds_rna)

# Gets normalization factors and normalizes counts in the sample comparison. Lastly,
# writes out normalized counts to a csv file.
# CLIP.
normalizationFactors <- normalizationFactors(dds)
clipNormalizedCounts<-counts(dds, normalized=TRUE)
write.csv(as.data.frame(clipNormalizedCounts), 
          file="all_clip_normalized_counts.csv")
# CLIP/SM.
normalizationFactors <- normalizationFactors(dds)
clipNormalizedCounts<-counts(dds, normalized=TRUE)
write.csv(as.data.frame(clipNormalizedCounts), 
          file="all_clip_sm_normalized_counts.csv")
# RNA-Seq
normalizationFactors <- normalizationFactors(dds)
clipNormalizedCounts<-counts(dds, normalized=TRUE)
write.csv(as.data.frame(clipNormalizedCounts), 
          file="all_rnaseq_normalized_counts.csv")

# Compares each CLIP sample to its respective SM sample. The events that exceed the
# cutoff (adjusted P-value of 0.10 and an absolute value of log2 fold change of 1.00)
# are considered true binding regions. These regions are used to filter in the
# normalized counts file and the concentration comparison files.
# 50x MBNL1 samples.
x50_vs_x50_sm <-results(dds_clip_sm, contrast=c("condition","x50","x50_sm"))
(x50_vs_x50_sm_ordered <- x50_vs_x50_sm [order(x50_vs_x50_sm$padj), ])
x50_vs_x50_sm_filtered <- subset(x50_vs_x50_sm_ordered, padj < 0.10)
x50_vs_x50_sm_filtered <- subset(x50_vs_x50_sm_filtered, log2FoldChange >= 1.00)
write.csv(as.data.frame(x50_vs_x50_sm_filtered), 
          file="x50_clip_sm_filtered.csv")
# 36x MBNL1 samples.
x36_vs_x36_sm <-results(dds_clip_sm, contrast=c("condition","x36","x36_sm"))
(x36_vs_x36_sm_ordered <- x36_vs_x36_sm [order(x36_vs_x36_sm$padj), ])
x36_vs_x36_sm_filtered <- subset(x36_vs_x36_sm_ordered, padj < 0.10)
x36_vs_x36_sm_filtered <- subset(x36_vs_x36_sm_filtered, log2FoldChange >= 1.00)
write.csv(as.data.frame(x36_vs_x36_sm_filtered), 
          file="x36_clip_sm_filtered.csv")
# 20x MBNL1 samples.
x20_vs_x20_sm <-results(dds_clip_sm, contrast=c("condition","x20","x20_sm"))
(x20_vs_x20_sm_ordered <- x20_vs_x20_sm [order(x20_vs_x20_sm$padj), ])
x20_vs_x20_sm_filtered <- subset(x20_vs_x20_sm_ordered, padj < 0.10)
x20_vs_x20_sm_filtered <- subset(x20_vs_x20_sm_filtered, log2FoldChange >= 1.00)
write.csv(as.data.frame(x20_vs_x20_sm_filtered), 
          file="x20_clip_sm_filtered.csv")
# 10x MBNL1 samples.
x10_vs_x10_sm <-results(dds_clip_sm, contrast=c("condition","x10","x10_sm"))
(x10_vs_x10_sm_ordered <- x10_vs_x10_sm [order(x10_vs_x10_sm$padj), ])
x10_vs_x10_sm_filtered <- subset(x10_vs_x10_sm_ordered, padj < 0.10)
x10_vs_x10_sm_filtered <- subset(x10_vs_x10_sm_filtered, log2FoldChange >= 1.00)
write.csv(as.data.frame(x10_vs_x10_sm_filtered), 
          file="x10_clip_sm_filtered.csv")
# 5x MBNL1 samples.
x5_vs_x5_sm <-results(dds_clip_sm, contrast=c("condition","x5","x5_sm"))
(x5_vs_x5_sm_ordered <- x5_vs_x5_sm [order(x5_vs_x5_sm$padj), ])
x5_vs_x5_sm_filtered <- subset(x5_vs_x5_sm_ordered, padj < 0.10)
x5_vs_x5_sm_filtered <- subset(x5_vs_x5_sm_filtered, log2FoldChange >= 1.00)
write.csv(as.data.frame(x5_vs_x5_sm_filtered), 
          file="x5_clip_sm_filtered.csv")
# Uninduced MBNL1 samples.
uni_vs_uni_sm <-results(dds_clip_sm, contrast=c("condition","uni","uni_sm"))
(uni_vs_uni_sm_ordered <- uni_vs_uni_sm [order(uni_vs_uni_sm$padj), ])
uni_vs_uni_sm_filtered <- subset(uni_vs_uni_sm_ordered, padj < 0.10)
uni_vs_uni_sm_filtered <- subset(uni_vs_uni_sm_filtered, log2FoldChange >= 1.00)
write.csv(as.data.frame(uni_vs_uni_sm_filtered), 
          file="uni_clip_sm_filtered.csv")

# Compares each induced MBNL1 sample set to the uninduced sample set. This works as
# an additional filtering mechanism. This filters for only binding regions that increase
# in binding as MBNL1 is induced. Filtering uses adjusted P-value of 0.10 and an 
# absolute value of log2 fold change of 0.5.
# Uninduced versus 50x.
x50_vs_uni <-results(dds_clip, contrast=c("condition","x50","uni"))
(x50_vs_uni_ordered <- x50_vs_uni[order(x50_vs_uni$padj), ])
x50_vs_uni_filtered <- subset(x50_vs_uni_ordered, padj < 0.10)
x50_vs_uni_filtered <- subset(x50_vs_uni_filtered, log2FoldChange > 0.5 | log2FoldChange < -0.5)
write.csv(as.data.frame(x50_vs_uni_filtered), 
          file="x50_vs_uni_significant.csv")
# Uninduced versus 36x.
x36_vs_uni <-results(dds_clip, contrast=c("condition","x36","uni"))
(x36_vs_uni_ordered <- x36_vs_uni[order(x36_vs_uni$padj), ])
x36_vs_uni_filtered <- subset(x36_vs_uni_ordered, padj < 0.10)
x36_vs_uni_filtered <- subset(x36_vs_uni_filtered, log2FoldChange > 0.5 | log2FoldChange < -0.5)
write.csv(as.data.frame(x36_vs_uni_filtered), 
          file="x36_vs_uni_significant.csv")
# Uninduced versus 20x.
x20_vs_uni <-results(dds_clip, contrast=c("condition","x20","uni"))
(x20_vs_uni_ordered <- x20_vs_uni[order(x20_vs_uni$padj), ])
x20_vs_uni_filtered <- subset(x20_vs_uni_ordered, padj < 0.10)
x20_vs_uni_filtered <- subset(x20_vs_uni_filtered, log2FoldChange > 0.5 | log2FoldChange < -0.5)
write.csv(as.data.frame(x20_vs_uni_filtered), 
          file="x20_vs_uni_significant.csv")
# Uninduced versus 10x.
x10_vs_uni <-results(dds_clip, contrast=c("condition","x10","uni"))
(x10_vs_uni_ordered <- x10_vs_uni[order(x10_vs_uni$padj), ])
x10_vs_uni_filtered <- subset(x10_vs_uni_ordered, padj < 0.10)
x10_vs_uni_filtered <- subset(x10_vs_uni_filtered, log2FoldChange > 0.5 | log2FoldChange < -0.5)
write.csv(as.data.frame(x10_vs_uni_filtered), 
          file="x10_vs_uni_significant.csv")
# Uninduced versus 5x.
x5_vs_uni <-results(dds_clip, contrast=c("condition","x5","uni"))
(x5_vs_uni_ordered <- x5_vs_uni[order(x5_vs_uni$padj), ])
x5_vs_uni_filtered <- subset(x5_vs_uni_ordered, padj < 0.10)
x5_vs_uni_filtered <- subset(x5_vs_uni_filtered, log2FoldChange > 0.5 | log2FoldChange < -0.5)
write.csv(as.data.frame(x5_vs_uni_filtered), 
          file="x5_vs_uni_significant.csv")

# Optional steps. Makes MA plot with all datapoints. Each protein concentration gets colored differently
# in plot.
# Makes general MA plot with the points from the highest protein concentration versus uninduced. This is
# because this comparison shows the most data points.
res <- lfcShrink(dds, contrast=c("condition","x50","uni")) 
with(res, plot(log2FoldChange, -log10(pvalue), pch=16, cex.lab = 1.25, cex.axis = 1.75, cex=0.6, main="Novel Regions Per Dose of MBNL1 Volcano Plot", cex.main = 2.5, 
               xlim=c(0,7)))
# Binding region name files of SM filtered protein concentration vs. uninduced. These files will have to be produced separately, after the
# normalized gene count files have been filtered with the SM comparison binding regions. The files contain one untitled column with
# the binding region. It is possible to use the filtered counts files themselves, but I found this way easier at the time.
regions_x5_vs_uni <- scan("x5_vs_uni_significant_sm_filtered.csv", what="", sep=",")
regions_x10_vs_uni <- scan("x10_vs_uni_significant_sm_filtered.csv", what="", sep=",")
regions_x20_vs_uni <- scan("x20_vs_uni_significant_sm_filtered.csv", what="", sep=",")
regions_x36_vs_uni <- scan("x36_vs_uni_significant_sm_filtered.csv", what="", sep=",")
regions_x50_vs_uni<- scan("x50_vs_uni_significant_sm_filtered.csv", what="", sep=",")
regions_x5_vs_uni_vector <- as.vector(regions_x5_vs_uni)
regions_x10_vs_uni_vector <- as.vector(regions_x10_vs_uni)
regions_x20_vs_uni_vector <- as.vector(regions_x20_vs_uni)
regions_x36_vs_uni_vector <- as.vector(regions_x36_vs_uni)
regions_x50_vs_uni_vector <- as.vector(regions_x50_vs_uni)

# Adds colors to MA plot. The different colors indicate the signifcant binding regions at the first protein
# concentration they become significant.
with(subset(res, rownames(res) %in% gene_list), points(log2FoldChange, -log10(pvalue), pch=10, col="blue"))
with(subset(res, rownames(res) %in% gene_list1), points(log2FoldChange, -log10(pvalue), pch=10, col="red"))
with(subset(res, rownames(res) %in% gene_list2), points(log2FoldChange, -log10(pvalue), pch=10, col="green"))
with(subset(res, rownames(res) %in% gene_list3), points(log2FoldChange, -log10(pvalue), pch=10, col="gold"))
with(subset(res, rownames(res) %in% gene_list4), points(log2FoldChange, -log10(pvalue), pch=10, col="purple"))

# Makes separate MA plots with signifcant binding regions highlighted in red.
# Uninduced versus 50x.
resLFC <- lfcShrink(dds, coef="condition_x50_vs_uni", type="apeglm")
plotMA(resLFC, ylim=c(-3,7), main="Uni vs 50x MBNL1",
       cex.main=2.25, cex.lab=1.5, cex.axis=1.5, xlab="Mean of Normalized Counts", ylab="Log2 Fold Change",) 
# Uninduced versus 36x.
resLFC <- lfcShrink(dds, coef="condition_x36_vs_uni", type="apeglm")
plotMA(resLFC, ylim=c(-3,7), main="Uni vs 36x MBNL1",
       cex.main=2.25, cex.lab=1.5, cex.axis=1.5, xlab="Mean of Normalized Counts", ylab="Log2 Fold Change",) 
# Uninduced versus 20x.
resLFC <- lfcShrink(dds, coef="condition_x20_vs_uni", type="apeglm")
plotMA(resLFC, ylim=c(-3,7), main="Uni vs 20x MBNL1",
       cex.main=2.25, cex.lab=1.5, cex.axis=1.5, xlab="Mean of Normalized Counts", ylab="Log2 Fold Change",) 
# Uninduced versus 10x.
resLFC <- lfcShrink(dds, coef="condition_x10_vs_uni", type="apeglm")
plotMA(resLFC, ylim=c(-3,7), main="Uni vs 10x MBNL1",
       cex.main=2.25, cex.lab=1.5, cex.axis=1.5, xlab="Mean of Normalized Counts", ylab="Log2 Fold Change",) 
# Uninduced versus 5x.
resLFC <- lfcShrink(dds, coef="condition_x5_vs_uni", type="apeglm")
plotMA(resLFC, ylim=c(-3,7), main="Uni vs 5x MBNL1",
       cex.main=2.25, cex.lab=1.5, cex.axis=1.5, xlab="Mean of Normalized Counts", ylab="Log2 Fold Change",) 



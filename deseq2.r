library("DESeq2")

counts <- as.matrix(read.csv(file = "WT_Cold3h_STRIPE-Seq_featureCounts_Ensembl.txt", sep="\t", row.names = "Geneid")) 
coldata <- read.csv(file = "condition.txt", sep="\t", row.names=1)
matrix <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design =~ condition)
                              
dds <- DESeq(matrix)

res <- results(dds, contrast=c("condition","cold","WT"), lfcThreshold=1.0, alpha=0.05 )
resSig <- subset(res, padj < 0.1)

summary(res)
summary(resSig)

write.csv(as.data.frame(res),file = "DESeq_STRIPE-Seq_cold3h_WT_genes_Ensembl.csv")
write.csv(as.data.frame(resSig),file = "DESeq_STRIPE-Seq_cold3h_WT_genes_padjLT0.1_Ensembl.csv")


de_res <- read.table("DESeq_STRIPE-Seq_cold3h_WT_genes_Ensembl.csv", sep = ",", dec = ".", header=TRUE)
head(de_res)

# Make a basic volcano plot
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.0, cex.lab=1.4, lwd=2.0)

with(de_res, plot(log2FoldChange, -log10(pvalue), pch=20, cex = 2.0, main="STRIPE-Seq WT vs Cold_3h Ensembl", xlim=c(-10,10)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(de_res, padj<.1 & log2FoldChange<=-1), points(log2FoldChange, -log10(pvalue), pch=20, cex = 2.0, col="blue"))
with(subset(de_res, padj<.1 & log2FoldChange>=1), points(log2FoldChange, -log10(pvalue), pch=20, cex = 2.0, col="red"))

abline(v=-1, col="black", lty=4, lwd=2.0)
abline(v=1, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(de_res$pvalue[de_res$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

dev.copy2pdf(file = "STRIPE-Seq_WTvsCold3h_Volcano_Up12_Down5_Ensembl.pdf")
 

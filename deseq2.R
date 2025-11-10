# script used for DESeq2 analysis
library(DESeq2)

count_data = as.matrix(read.csv("MF_gonads_FCforDS.tsv",sep="\t",row.names="gene_id"))
count_data <- count_data[, !colnames(count_data) %in% "FOS84"] #P24XY samples removed after PCA analysis
col_data = read.table(file = "MF_gonads_info.txt", header = T, sep = "\t")

count_data = as.matrix(read.csv("MF_head_FCforDS.tsv",sep="\t",row.names="gene_id"))
count_data <- count_data[, !colnames(count_data) %in% "FHS8"] #P24X0 samples removed after PCA analysis
count_data <- count_data[, !colnames(count_data) %in% "FHS61"] #P24XY samples removed after PCA analysis
col_data = read.table(file = "MF_head_info.txt", header = T, sep = "\t")

count_data = as.matrix(read.csv("MF_legs_FCforDS.tsv",sep="\t",row.names="gene_id"))
remove_samples <- c("FLS20", "MLS39") #P24X0 samples removed after PCA analysis
count_data <- count_data[, !colnames(count_data) %in% remove_samples] 
count_data <- count_data[, !colnames(count_data) %in% "FLS74"] #P24XY samples removed after PCA analysis
col_data = read.table(file = "MF_legs_info.txt", header = T, sep = "\t")

dds = DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "FOS") # changed to FHS for head samples & FLS for leg samples respectively
dds = DESeq(dds)
bmab <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) ) # to get base expression of male and female samples
res <- results(dds, alpha=0.05, lfcThreshold=2) # padj and log2FC cutoffs
z <- cbind(gene_id = rownames(res), bmab, res)
z$padj[!is.na(z$pvalue) & z$pvalue == 1] <- 1 # to ensure genes are not treated as false-positives
z_filtered_male <- subset(z, padj < 0.05 & log2FoldChange > 2) # male-biased genes
z_filtered_female <- subset(z, padj < 0.05 & log2FoldChange < -2) # female-biased genes

# for P24X0 race, changed gonads to head / legs respectively
write.table(z, file="MF_gonads_onP24XOM_DEGs.tsv", sep='\t', quote=F, row.names=F)
write.table(z_filtered_male, file="MF_gonads_onP24XOM_mbDEGs.tsv", sep='\t', quote=F, row.names=F)
write.table(z_filtered_female, file="MF_gonads_onP24XOM_fbDEGs.tsv", sep='\t', quote=F, row.names=F)

# for P24XY race, changed gonads to head / legs respectively
write.table(z, file="MF_gonads_onP24XYF_DEGs.tsv", sep='\t', quote=F, row.names=F)
write.table(z_filtered_male, file="MF_gonads_onP24XYF_mbDEGs.tsv", sep='\t', quote=F, row.names=F)
write.table(z_filtered_female, file="MF_gonads_onP24XYF_fbDEGs.tsv", sep='\t', quote=F, row.names=F)

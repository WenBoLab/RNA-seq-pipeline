# Manuscript:   Downstream analysis of RNA-seq data: 
#               Differential expression analysis、Enrichment analysis and WGCNA
# Author:       Shaolizhen
#=======================================================================================

# load libraries
library(DESeq2)   
library(ggplot2)
library(tidyverse)
library(tidyr)
library(ggrepel)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)
library(tibble)
library(WGCNA)
library(VennDiagram)
library(org.Mm.eg.db)
library(ReactomePA)
library(scales)
library(stringr)
library(msigdbr)
library(enrichplot)
# set working directory
data_dir <- "E:/RNA-seq/01/data"
result_dir <- 'E:/RNA-seq/01/results'
# set working directory to current file location
setwd(data_dir)
# step1 make source data -------------------------------------------------------------------
# input data
countData <- read.table("gene_counts.txt",sep = "\t", header = T)
# check data
head(countData)
# change sample names
sampleNames = c("RNA12","RNA13","RNA14","RNA15","RNA16","RNA17","RNA18","RNA19")
names(countData)[7:14] = sampleNames
# remove low expression genes
countData <- countData[rowSums(countData[c(7:14)]) > 0, ]   # 38956
# convert geneid to gene name
gene_info <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene2symbols <- getBM(attributes = c('external_gene_name','ensembl_gene_id','gene_biotype'),
                      filters = "ensembl_gene_id", values = countData$Geneid, mart = gene_info)
unique(gene2symbols$gene_biotype)

countData.id = merge(gene2symbols, countData,by.x='ensembl_gene_id',by.y='Geneid')

# Set gene names as rownames
rownames(countData.id) <- countData.id$external_gene_name
# remove duplicated genes
is_duplicate <- duplicated(countData.id$external_gene_name)
data <- countData.id[!is_duplicate, ]
rownames(data) <- data$external_gene_name
data <- data[,-c(1:5)]

# TPM data
gene_length_kb <- data$Length / 1000
count_mat <- as.matrix(data[ , 2:ncol(data)])
rpk <- sweep(count_mat, 1, gene_length_kb, FUN = "/")
sum_rpk <- colSums(rpk)
tpm <- as.data.frame(sweep(rpk, 2, sum_rpk, FUN = "/") * 1e6)

# save data
setwd(result_dir)
write.table(tpm,file = "1.data_with_genes_tpm.txt",row.names = TRUE,col.names = TRUE,sep = '\t')

data = data[,-1]
write.table(data,file = "1.data_with_genes_counts.txt",row.names = TRUE,col.names = TRUE,sep = '\t')

# step2 differential expression analysis -------------------------------------------------------------------
# set parameters
cut_off_pvalue = 0.05
cut_off_logFC = 0.5
# set working directory to save file
setwd(result_dir)
# group
group = factor(c(rep("control",4), rep("treat",4)))
colData <- data.frame(group)
colData <- data.frame(row.names=colnames(data), group)

dds <- DESeqDataSetFromMatrix(countData = data, 
                              colData = colData, 
                              design = ~group)
dds <- DESeq(dds)

# save data
save(dds, file = '2.diff_RNA_dds.RData')

diff_data <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table
# sort with adjust pvalue
diff_data <- diff_mus_liver[order(diff_data$padj),]
head(diff_data)
diff_data <- as.data.frame(diff_data)
diff_data$gene <- rownames(diff_data)

write.table(diff_data,file = "2.diff_genes.txt",row.names = TRUE,col.names = TRUE,sep = '\t')

up_genes <- diff_data[which(diff_data$pvalue < cut_off_pvalue & diff_data$log2FoldChange > cut_off_logFC),]
down_genes <- diff_data[which(diff_data$pvalue < cut_off_pvalue & diff_data$log2FoldChange < -cut_off_logFC),]
write.table(up_genes,file = "2.diff_genes_UP.txt",row.names = TRUE,col.names = TRUE,sep = '\t')
write.table(down_genes,file = "2.diff_genes_DOWN.txt",row.names = TRUE,col.names = TRUE,sep = '\t')

# step3 correlation between samples -------------------------------------------------------------------
rld <- rlog(dds, blind=TRUE)  # For small sample sizes, the rlog method is recommended; for large sample sizes, the vst method is preferred.
mat_rld <- assay(rld)

pdf('2.RNA-sample_hclust.pdf',width = 4,height = 4)
dists <- dist(t(assay(rld)))
plot(hclust(dists, method="ward.D2"))
dev.off()

pdf('2.RNA-PCA.pdf',width = 6,height = 6)
p <- plotPCA(rld, intgroup="group")
p + geom_text_repel(aes(label = name), size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("control" = "#4575B4", "treat" = "#D73027")) +
  theme_bw()
dev.off()

cor_mat <- cor(mat_rld)
pdf('2.RNA-correlation.pdf',width = 6,height = 5)
pheatmap(cor_mat,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Sample correlation (rlog transformed)")
dev.off()

# step4 cvolcanoplot of diff genes ------------------------------------------------------------------- 
apply(diff_data,2, function(x) sum(is.na(x)))  # count NA
diff_data = drop_na(diff_data, pvalue)  # remove NA
apply(diff_data,2, function(x) sum(is.na(x)))  # again

# set gene status
diff_data$change = ifelse(diff_data$pvalue < cut_off_pvalue & abs(diff_data$log2FoldChange) >= cut_off_logFC, 
                               ifelse(diff_data$log2FoldChange > cut_off_logFC ,'Up','Down'),
                               'Stable')
table(diff_data$change)
# add gene labels in plot
data <- diff_data %>%
  dplyr::filter(padj < 0.05) %>% 
  arrange(log2FoldChange) 
data$label = ""
data$label[1:10] = data$gene[1:10]

data <- data %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(-log2FoldChange) 
data$label[1:10] = data$gene[1:10]

diff_data$label <- data$label[match(diff_data$gene, data$gene)]

p <- ggplot(
  diff_data, aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
  geom_point(alpha=0.8, size=1) +
  scale_color_manual(values=c("#104F8D", "#d2dae2","#BD4F48"))+
  geom_vline(xintercept=c(-log2FoldChange,log2FoldChange),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(Fold Change)",
       y="-log10 (pvalue)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) +
  geom_text_repel(aes(x = log2FoldChange,                        
                      y = -1*log10(pvalue),          
                      label=label),                       
                  max.overlaps = 10000,                     
                  size=2.5,                                  
                  box.padding=unit(0.5,'lines'),           
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',                    
                  show.legend=FALSE)

pdf("3.volcanoplot of diff genes.pdf", width = 6, height = 6)
p
dev.off()

# step5 heatmap of diff genes ------------------------------------------------------------------- 
diff_data_s <- diff_data[which(diff_data$change %in% c('Up','Down')),]
# use tpm to display
data_tpm_s <- tpm[which(rownames(tpm) %in% rownames(diff_data_s)),]

annotation_col <- rep(c(rep("control",4), rep("treat",4)))
annotation_col <- as.data.frame(annotation_col,row.names = colnames(data_tpm_s))

col_colors = c('control' = '#00DAE0', 'treat' = '#FF9289')

pdf('3.heatmap of diff genes.pdf',width = 6,height = 6)
pheatmap(data_tpm_s,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "row",
         legend = T,
         show_colnames = TRUE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = list(annotation_col = col_colors))
dev.off()

# step6 enrichment analysis of diff genes ------------------------------------------------------------------- 
pvalueSET = 0.05 
qvalueSET = 0.05 
showNum = 10 

# UP
gene_info <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene2symbols <- getBM(attributes = c('external_gene_name','entrezgene_id'),
                      filters = "external_gene_name", values = UP_genes$gene, mart = gene_info)

go <- enrichGO(gene = gene2symbols$entrezgene_id,
               OrgDb = org.Mm.eg.db,     
               keyType = "ENTREZID",     
               ont = "ALL",              
               pAdjustMethod = "BH",     
               pvalueCutoff = pvalueSET, 
               qvalueCutoff = qvalueSET, 
               readable = T             
)
go.data <- data.frame(go) 
write.table(go.data,file="3.GO of UP genens.txt",sep="\t",quote=F,row.names = TRUE)

go.data.bp <- go.data[go.data$`ONTOLOGY` == "BP", ]
go.data.cc <- go.data[go.data$`ONTOLOGY` == "CC", ]
go.data.mf <- go.data[go.data$`ONTOLOGY` == "MF", ]

go.data.bp <- go.data.bp[order(go.data.bp$`pvalue`), ]
go.data.cc <- go.data.cc[order(go.data.cc$`pvalue`), ]
go.data.mf <- go.data.mf[order(go.data.mf$`pvalue`), ]

go.data.bp.top <- head(go.data.bp, n = showNum)
go.data.cc.top <- head(go.data.cc, n = showNum)
go.data.mf.top <- head(go.data.mf, n = showNum)


go.data.top <- rbind(go.data.bp.top, go.data.cc.top, go.data.mf.top)
go.data.top$Description <- factor(go.data.top$Description,levels = rev(go.data.top$Description))

gene_ratio_num <- sapply(strsplit(as.character(go.data.top$GeneRatio), "/"),
                         function(x) as.numeric(x[1]) / as.numeric(x[2]))

go.data.top$GeneRatio_num <- gene_ratio_num

p <- ggplot(go.data.top, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  coord_cartesian(clip = "off") +
  scale_color_paletteer_c(
    palette   = "viridis::plasma",
    direction = -1,
    name      = "p.adjust"
  ) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 5) 
  ) +
  labs(
    x = "GeneRatio",
    y = "Terms"
  ) +
  facet_grid(ONTOLOGY ~., scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 9)
  )
pdf("3.GO dotplot of UP genes.pdf", width = 8, height = 8)
p
dev.off()

kegg_enrich <- enrichKEGG(
  gene = gene2symbols$entrezgene_id,
  keyType = 'kegg',
  organism = 'mmu',
  pvalueCutoff = pvalueSET,
  qvalueCutoff = qvalueSET
)

kegg.data <- setReadable(kegg_enrich, 'org.Mm.eg.db', 'ENTREZID')
kegg.data <- as.data.frame(kegg.data)

#kegg.data$Description <- str_split(kegg.data$Description,' - ',simplify = T)[,1]
write.table(kegg.data,file="3.KEGG of UP genens.txt",sep="\t",quote=F,row.names = TRUE)
rownames(kegg.data) <- 1:nrow(kegg.data)
kegg.data$order <- factor(rev(as.integer(rownames(kegg.data))),labels = rev(kegg.data$Description))

kegg_bar <- ggplot(data = kegg.data, aes(y = order, x = Count)) +
  geom_point(aes(size=Count,color=p.adjust)) + 
  scale_color_paletteer_c(
    palette   = "viridis::plasma",
    direction = -1,
    name      = "p.adjust"
  ) +
  facet_grid(category ~., scales = "free_y",space = "free_y") +
  theme_bw() + # 黑白主题
  labs(x = "Gene Counts", y = "Pathways", title = "KEGG Enrichment") +
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"), 
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

pdf("3.KEGG dotplot of UP genes.pdf", width = 6, height = 5)
kegg_bar
dev.off()


# the downregulated genes enrichment analysis as same as upregulated genes


# step7 GSEA analysis ------------------------------------------------------------------- 

lfc_vector <- diff_data$log2FoldChange
names(lfc_vector) <- diff_data$gene

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
head(lfc_vector)

geneset <- msigdbr(species = "Mus musculus",category = "H")

# start GSEA
set.seed(123)
gsea_mm <- GSEA(geneList = lfc_vector,
                minGSSize = 10,                  
                maxGSSize = 500,                   
                pvalueCutoff = 0.25,                    
                pAdjustMethod= "BH",                  
                seed = TRUE,
                TERM2GENE = dplyr::select(
                  geneset,
                  gs_name,
                  gene_symbol
                ))  

head(gsea_mm@result)
save(gsea_mm, file = '4.GSEA_ALLgenes.RData')

gsea_result_df <- data.frame(gsea_mm@result)
write.table(gsea_result_df, "4.GSEA_Hallmark_ALL.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_mm,
  geneSetID = "HALLMARK_PROTEIN_SECRETION",
  title = "HALLMARK_PROTEIN_SECRETION",
  color.line = "#0d76ff"
)
most_positive_nes_plot

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

most_negative_nes_plot <- enrichplot::gseaplot2(
  gsea_mm,
  geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  title = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
most_negative_nes_plot


pdf('GSEA_Hallmark_split.pdf', width = 6, height = 6)
lapply(1:nrow(gsea_mm@result), function(x) {
  gseaplot2(gsea_mm, geneSetID = gsea_mm[x]$Description[1], title = gsea_mm[x]$Description[1])
})
dev.off()

pdf('GSEA_Hallmark_ALL.pdf', width = 5, height = 4)
dotplot(gsea_mm, showCategory = 8, orderBy = "p.adjust") +
  scale_color_gradient(low = "red", high = "blue") +
  scale_size(range = c(3, 8)) +
  theme_bw()

dev.off()

# step8 WGCNA analysis ------------------------------------------------------------------- 
# referencr https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#47_Run_WGCNA
#           https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
networkType <- "unsigned"  # "unsigned" 或 "signed"-方向一致

sample_info <- data.frame(group)
sample_info <- data.frame(row.names=colnames(data), group)

data <- data[rowSums(data) > 10, ]
# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(countData = data, 
                              colData = sample_info, 
                              design = ~1  # Here we are not specifying a model
)
dds_norm <- vst(dds)
plotPCA(dds_norm, intgroup="group")
# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>% t()

gsg <- goodSamplesGenes(normalized_counts, verbose = 3)
gsg$allOK
if (!gsg$allOK) {   
  normalized_counts <- normalized_counts[, gsg$goodGenes]
  normalized_counts <- normalized_counts[gsg$goodSamples, ]
}
summary(gsg)

sampleTree <- hclust(dist(normalized_counts), method = "average")
plot(sampleTree, main="Sample clustering", sub="", xlab="")

sft <- pickSoftThreshold(normalized_counts, 
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = networkType
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

plot <- ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

pdf(file = '5.WGCNA_pickSoftThreshold.pdf', width = 4,height = 4)
print(plot)
dev.off()

#softPower <- sft$powerEstimate
softPower <- 20
message(sprintf("Chosen soft-thresholding power = %s", softPower))

bwnet <- WGCNA::blockwiseModules(
  normalized_counts,
  power = softPower,
  networkType = networkType,  # "unsigned" / "signed"
  TOMType = networkType,
  maxBlockSize = 10000,        
  numericLabels = TRUE,
  randomSeed = 1234
)

save(bwnet, file = '5.WGCNA_TOM.RData')

module_eigengenes <- bwnet$MEs
head(module_eigengenes)

gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))
write.table(gene_module_key,file = "5.WGCNA_gene_module_key.txt",row.names = TRUE,col.names = TRUE,sep = '\t')
write.table(module_eigengenes,file = "5.WGCNA_module_eigengenes.txt",row.names = TRUE,col.names = TRUE,sep = '\t')


all.equal(rownames(sample_info), rownames(module_eigengenes))

pdf(file = '5.WGCNA_group.pdf', width = 4,height = 3)
des_mat <- model.matrix(~ sample_info[,1])
# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)
# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")
write.table(stats_df,file = "4.WGCNA_module_eigengenes_stats_df.txt",row.names = TRUE,col.names = TRUE,sep = '\t')

head(stats_df)
sample_info$sample_id <- rownames(sample_info)
modules <- stats_df[which(stats_df$P.Value<0.05),1] 
module_eigengenes <- module_eigengenes[,which(colnames(module_eigengenes) %in% modules)]

module_df <- module_eigengenes %>%
  tibble::rownames_to_column("accession_code") %>%
  # Here we are performing an inner join with a subset of sample_info
  dplyr::inner_join(sample_info %>%
                      dplyr::select(sample_id, group),
                    by = c("accession_code" = "sample_id")
  )  
for (j in 1:length(module_eigengenes)){
  module = names(module_eigengenes[j])
  p <- ggplot(module_df,aes(x = sample_info[,1],y = module_eigengenes[,j],color = sample_info[,1])) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggforce::geom_sina(maxwidth = 0.3) +
    theme_classic() +
    ggtitle(module) +
    xlab('Group') +
    ylab(module) +
    labs(color = "Group") 
  print(p)
}

dev.off()

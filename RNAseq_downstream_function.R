#' @description   RNA-seq pipeline
#' @title  1）差异分析
#'   2）富集分析

#' @author  shaolizhen
Plus.library <- function(packages){
  ## 系统已经安装的包
  x <- installed.packages()
  x <- as.data.frame(x)
  installed <- as.character(x$Package)
  
  ## 检查系统中有没有安装pacman包
  is <- length(grep("pacman",installed))
  if(is ==0){
    #pacman包未安装
    install.packages("pacman")
    library(pacman)
  } else {
    library(pacman)
  }
  
  ## 加载包
  for(i in packages){p_load(char=i)}
}

#' @description 判断基因编码 ID 类型
#' @param gene_ids 基因编码 ID 
#' @return 基因编码类型
detect_gene_id_type <- function(gene_ids, prop_threshold = 0.9) {
  gene_ids <- gene_ids[!is.na(gene_ids)]
  is_ens_version  <- grepl("^ENS[G|T][0-9]+\\.[0-9]+$", gene_ids)
  is_ens_stable   <- grepl("^ENS[G|T][0-9]+$", gene_ids) & !grepl("\\.", gene_ids)
  is_symbol_like  <- !is_ens_version & !is_ens_stable
  
  prop_ens_version <- mean(is_ens_version)
  prop_ens_stable  <- mean(is_ens_stable)
  prop_symbol      <- mean(is_symbol_like)
  
  # 调试用（可注释掉）
  # message(sprintf("比例: ens.version=%.2f%%, ens.stable=%.2f%%, symbol/other=%.2f%%", 
  #                 prop_ens_version*100, prop_ens_stable*100, prop_symbol*100))
  
  if (prop_ens_version >= prop_threshold) {
    return("ensembl_gene_id_version")
  }
  if (prop_ens_stable >= prop_threshold) {
    return("ensembl_gene_id")
  }
  if (prop_symbol >= prop_threshold) {
    return("external_gene_name")
  }
  
  # 实在混杂得很严重
  return("mixed_or_unrecognized")
}

#' @description 基因 ID 转化
#' @param gene_ids 基因编码 ID 
#' @return 基因信息
convert_to_gene_symbol <- function(gene_ids,
                                   species = c("human", "mouse")) {
  # import library
  if(species == "human"){
    nd <- c("biomaRt","org.Hs.eg.db")
    Plus.library(nd)
    bio.use = "hsapiens_gene_ensembl"
  }
  
  if(species == "mouse"){
    nd <- c("biomaRt","org.Mm.eg.db")
    Plus.library(nd)
    bio.use = "mmusculus_gene_ensembl"
  }
  id_type <- detect_gene_id_type(gene_ids)
  # gene2symbol
  gene_info <- useDataset(bio.use, useMart("ensembl"))
  gene2symbols <- getBM(attributes = c('external_gene_name','entrezgene_id','description','gene_biotype',id_type),
                        filters = id_type, values = gene_ids, mart = gene_info)
  
  return(gene2symbols)
}

check_is_raw_counts <- function(count_mat) {
  if (!is.matrix(count_mat) && !is.data.frame(count_mat)) {
    stop("Input count data must be a matrix or data.frame.")
  }
  
  if (any(count_mat < 0, na.rm = TRUE)) {
    stop("Count matrix contains negative values.")
  }
  
  if (!all(count_mat == round(count_mat), na.rm = TRUE)) {
    warning(
      "Input data does NOT appear to be raw counts (non-integer values detected).\n",
      "DESeq2 requires raw integer counts. Please check your input."
    )
    return(FALSE)
  }
  message("✔ Raw count data detected.")
  return(TRUE)
}

#' @description 差异分析
#' @param data.dir 定量数据存放地址
#' @param count_file raw count文件名
#' @param sample_info_file 样本信息
run_DESeq2 <- function(
    data.dir,
    count_file,
    sample_info_file,
    group_col = "condition",
    control_level,
    species = c("human", "mouse"),
    filter_low_counts = TRUE,
    min_counts = 0,
    min_samples = 0,
    select = c('padj','pvalue'),
    cut_off_pvalue = 0.05,
    cut_off_logFC = 1) {
  
  # import library
  nd <- c("DESeq2","dplyr","readxl",'tidyr')
  Plus.library(nd)
  # read data
  counts <- read.table(file=file.path(data.dir,count_file), header = TRUE, row.names = 1)
  # filter data
  if(filter_low_counts == TRUE){
    counts <- counts[rowMeans(counts[,6:ncol(counts)]) > min_counts, ]
  }
  # sample information
  sample_info <- read_excel(file.path(data.dir,sample_info_file))
  sample_info <- as.data.frame(sample_info)
  rownames(sample_info) <- sample_info$sample_name
  
  ## rename sample
  new_names <- paste0(sample_info$condition,"_",sample_info$replicate)
  colnames(counts)[6:ncol(counts)] <- new_names
  check_is_raw_counts(counts[,6:ncol(counts)])
  
  ## gene annotation
  gene_id_type <- detect_gene_id_type(rownames(counts))
  gene_anno <- convert_to_gene_symbol(rownames(counts), species)
  
  ## prepare metadata
  sample_info[[group_col]] <- factor(sample_info[[group_col]])
  sample_info[[group_col]] <- relevel(sample_info[[group_col]], ref = control_level)
  ## DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = counts[, 6:ncol(counts)],
    colData = sample_info,
    design = as.formula(paste0("~", group_col))
  )

  dds <- DESeq(dds)
  res <- results(dds)
  
  if(select == 'pvalue'){
    diff_data <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(gene_anno, by = setNames(gene_id_type, "gene_id")) %>%
      mutate(
        change = case_when(
          pvalue < cut_off_pvalue & log2FoldChange >= cut_off_logFC ~ "Up",
          pvalue < cut_off_pvalue & log2FoldChange <= -cut_off_logFC ~ "Down",
          TRUE ~ "Stable"
        )
      ) %>%
      mutate(across(where(is.character), ~na_if(., ""))) %>%
      mutate(across(where(is.character), ~na_if(., "."))) 
    # remove NA
    if (any(is.na(diff_data$pvalue))) {
      cat("检测到 pvalue 中有 NA，正在移除...\n")
      diff_data <- drop_na(diff_data, pvalue)
    }
  }else{
    diff_data <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(gene_anno, by = setNames(gene_id_type, "gene_id")) %>%
      mutate(
        change = case_when(
          padj < cut_off_pvalue & log2FoldChange >= cut_off_logFC ~ "Up",
          padj < cut_off_pvalue & log2FoldChange <= -cut_off_logFC ~ "Down",
          TRUE ~ "Stable"
        )
      ) %>%
      mutate(across(where(is.character), ~na_if(., ""))) %>%
      mutate(across(where(is.character), ~na_if(., ".")))
    # remove NA
    if (any(is.na(diff_data$padj))) {
      cat("检测到 padj 中有 NA，正在移除...\n")
      diff_data <- drop_na(diff_data, padj)
    }
  }
  # remove NA in gene_name
  if (any(is.na(diff_data$external_gene_name))) {
    cat("检测到 gene_name 中有 NA，正在移除...\n")
    diff_data <- drop_na(diff_data, external_gene_name)
  }
  # remove duplicated gene name
  diff_data <- diff_data[!duplicated(diff_data$external_gene_name), ]
  rownames(diff_data) <- diff_data$external_gene_name
  
  return(list(dds = dds,DEG = diff_data))
}
plot_PCA <- function(dds,data.dir,sample_info_file){
  nd <- c("ggplot2","ggrepel")
  Plus.library(nd)
  sample_info <- read_excel(file.path(data.dir,sample_info_file))
  sample_info <- as.data.frame(sample_info)
  rownames(sample_info) <- sample_info$sample_name
  if(length(sample_info$sample_id) < 30){
    scale.data <- rlog(dds, blind=TRUE)
  }else{
    scale.data <- vst(dds, blind=TRUE)
  }
  p <- plotPCA(scale.data, intgroup="condition")
  return(p)
}

plot_volcanoplot <- function(diff_data,
                             cut_off_pvalue = 0.05,
                             cut_off_logFC = 1,
                             num = 5,
                             select = c('padj','pvalue')){
  nd <- c("ggplot2","paletteer","tidyverse","dplyr")
  select <- match.arg(select)
  cat('distribution of genes: \n')
  print(table(diff_data$change))
  # add gene labels in plot
  # ── 准备要标注的基因 ────────────────────────────────
  # 只考虑显著基因
  sig_data <- diff_data %>% 
    filter(padj < cut_off_pvalue)   # 或用 .data[[select]] < cut_off_pvalue
  
  # 上调（正 logFC）的前 num 个
  up <- sig_data %>%
    filter(log2FoldChange > cut_off_logFC) %>%
    arrange(desc(log2FoldChange)) %>%          # 从最大 logFC 开始
    slice_head(n = num) %>%
    mutate(direction = "Up")
  
  # 下调（负 logFC）的前 num 个
  down <- sig_data %>%
    filter(log2FoldChange < -cut_off_logFC) %>%
    arrange(log2FoldChange) %>%                # 从最小 logFC 开始
    slice_head(n = num) %>%
    mutate(direction = "Down")
  
  # 合并要标注的基因（最多 2*num 个）
  to_label <- bind_rows(up, down) %>%
    mutate(label = external_gene_name)
  
  # 把 label 对应回原数据框
  diff_data <- diff_data %>%
    left_join(
      to_label %>% select(external_gene_name, label),
      by = "external_gene_name"
    ) %>%
    mutate(
      label = replace_na(label, "")
    )
  
  p <- ggplot(
    diff_data, aes(x = log2FoldChange, y = -log10(.data[[select]]), colour=change)) +
    geom_point(alpha=0.8, size=1) +
    scale_color_manual(values=c("#104F8D", "#d2dae2","#BD4F48"))+
    geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
    labs(x="log2(Fold Change)",
         y=paste0("-log10 (",select,")"))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank()) +
    geom_text_repel(aes(label=label),                       
                    max.overlaps = 10000,                     
                    size=2.5,                                  
                    box.padding=unit(0.5,'lines'),           
                    point.padding=unit(0.1, 'lines'), 
                    segment.color='black',                    
                    show.legend=FALSE)
  
  return(p)
}

plot_heatmap <- function(
    data.dir,
    count_file,
    sample_info_file,
    gene_list,
    species = c("human", "mouse"),
    filter_low_counts = TRUE,
    min_counts = 0,
    min_samples = 0
) {
  nd <- c("pheatmap")
  Plus.library(nd)
  species <- match.arg(species)
  # read data
  counts <- read.table(file=file.path(data.dir,count_file), header = TRUE, row.names = 1)
  # filter data
  if(filter_low_counts == TRUE){
    counts <- counts[rowMeans(counts[,6:ncol(counts)]) > min_counts, ]
  }
  # sample information
  sample_info <- read_excel(file.path(data.dir,sample_info_file))
  sample_info <- as.data.frame(sample_info)
  rownames(sample_info) <- sample_info$sample_name
  
  ## rename sample
  new_names <- paste0(sample_info$condition,"_",sample_info$replicate)
  colnames(counts)[6:ncol(counts)] <- new_names
  
  ## gene annotation
  gene_id_type <- detect_gene_id_type(rownames(counts))
  gene_anno <- convert_to_gene_symbol(rownames(counts), species)
  counts <- counts %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(gene_anno, by = setNames(gene_id_type, "gene_id")) %>%
    mutate(across(where(is.character), ~na_if(., ""))) %>%
    mutate(across(where(is.character), ~na_if(., ".")))
  
  if (any(is.na(counts$external_gene_name))) {
    cat("检测到 gene_name 中有 NA，正在移除...\n")
    counts <- drop_na(counts, external_gene_name)
  }
  # remove duplicated gene name
  counts <- counts[!duplicated(counts$external_gene_name), ]
  rownames(counts) <- counts$external_gene_name
  
  if (is.null(counts$Length)) {
    counts <- counts[,rownames(sample_info)]
    counts <- counts[rownames(counts) %in% gene_list, ]
    apply(counts,2, function(x) sum(is.na(x)))  # count NA
    plot <- pheatmap(
      counts,
      main = "scale expression",
      scale = 'row',
      #annotation_col = sample_info$condition,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE
    )
  }
  else{
    gene_length_kb <- counts$Length / 1000
    rpk <- sweep(counts[, rownames(sample_info)], 1, gene_length_kb, FUN = "/")
    tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6
    tpm <- as.data.frame(tpm)
    
    tpm_sub <- tpm[rownames(tpm) %in% gene_list, ]
    
    plot <- pheatmap(
      tpm_sub,
      main = "TPM",
      scale = 'row,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE
    )
  }
  return(plot)
}

GO_enrichment <- function(gene_list,
                          species = c("human", "mouse"),
                          select = c('padj','pvalue'),
                          cut_off_pvalue = 0.05,
                          showNum = 10){
  
  pvalueSET = 1
  qvalueSET = 1
  ## gene annotation
  gene_id_type <- detect_gene_id_type(gene_list)
  gene_anno <- convert_to_gene_symbol(gene_list, species)
  if(species == "human"){
    nd <- c("clusterProfiler","org.Hs.eg.db","ggplot2","paletteer")
    Plus.library(nd)
    go <- enrichGO(gene = gene_anno$entrezgene_id,
                   OrgDb = org.Hs.eg.db,     
                   keyType = "ENTREZID",     
                   ont = "ALL",              
                   pAdjustMethod = "BH",     
                   pvalueCutoff = pvalueSET, 
                   qvalueCutoff = qvalueSET, 
                   readable = T             
    )
  }
  if(species == "mouse"){
    nd <- c("clusterProfiler","org.Mm.eg.db","ggplot2","paletteer")
    Plus.library(nd)
    go <- enrichGO(gene = gene_anno$entrezgene_id,
                   OrgDb = org.Mm.eg.db,     
                   keyType = "ENTREZID",     
                   ont = "ALL",              
                   pAdjustMethod = "BH",     
                   pvalueCutoff = pvalueSET, 
                   qvalueCutoff = qvalueSET, 
                   readable = T             
    )
  }
  
  go.data <- data.frame(go)
  if(select == 'pvalue'){
    go.data <- go.data[which(go.data$pvalue < cut_off_pvalue),]
    display_p <- 'pvalue'
  }else{
    go.data <- go.data[which(go.data$p.adjust < cut_off_pvalue),]
    display_p <- 'p.adjust'
  }
  
  # 根据生物过程（BP）、细胞组分（CC）和分子功能（MF）,分成三个表
  go.data.bp <- go.data[go.data$`ONTOLOGY` == "BP", ]
  go.data.cc <- go.data[go.data$`ONTOLOGY` == "CC", ]
  go.data.mf <- go.data[go.data$`ONTOLOGY` == "MF", ]
  
  # 按照p-value值升序对每个表格进行排序
  go.data.bp <- go.data.bp[order(go.data.bp$`pvalue`), ]
  go.data.cc <- go.data.cc[order(go.data.cc$`pvalue`), ]
  go.data.mf <- go.data.mf[order(go.data.mf$`pvalue`), ]
  # 从排序后的数据框中选择前20行 
  go.data.bp.top <- head(go.data.bp, n = showNum)
  go.data.cc.top <- head(go.data.cc, n = showNum)
  go.data.mf.top <- head(go.data.mf, n = showNum)
  
  # 将三个数据框合并成一个数据框
  go.data.top <- rbind(go.data.bp.top, go.data.cc.top, go.data.mf.top)
  go.data.top$Description <- factor(go.data.top$Description,levels = rev(go.data.top$Description))
  
  gene_ratio_num <- sapply(strsplit(as.character(go.data.top$GeneRatio), "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]))
  
  go.data.top$GeneRatio_num <- gene_ratio_num
  
  p <- ggplot(go.data.top, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = display_p)) +
    coord_cartesian(clip = "off") +
    # 颜色：用 viridis::plasma，并把颜色方向反过来（p越小颜色越亮）
    scale_color_paletteer_c(
      palette   = "viridis::plasma",
      direction = -1,              # 让小 p 在高亮端
      name      = display_p
    ) +
    scale_size_continuous(
      name = "Count",
      range = c(1, 5)             # 气泡稍微大一点更明显
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
  results <- list(data = go.data, plot=p)
  return(results)
}

KEGG_enrichment <- function(gene_list,
                          species = c("human", "mouse"),
                          select = 'padj',
                          cut_off_pvalue = 0.05,
                          showNum = 10){
  
  pvalueSET = 1
  qvalueSET = 1
  ## gene annotation
  gene_id_type <- detect_gene_id_type(gene_list)
  gene_anno <- convert_to_gene_symbol(gene_list, species)
  if(species == "human"){
    nd <- c("clusterProfiler","org.Hs.eg.db","ggplot2","paletteer")
    Plus.library(nd)
    kegg_enrich <- enrichKEGG(
      gene = gene_anno$entrezgene_id,
      keyType = 'kegg',
      organism = 'hsa',
      pvalueCutoff = pvalueSET,
      qvalueCutoff = qvalueSET
    )
    kegg.data <- setReadable(kegg_enrich, 'org.Hs.eg.db', 'ENTREZID')
  }
  if(species == "mouse"){
    nd <- c("clusterProfiler","org.Mm.eg.db","ggplot2","paletteer")
    Plus.library(nd)
    kegg_enrich <- enrichKEGG(
      gene = gene2symbols$entrezgene_id,
      keyType = 'kegg',
      organism = 'mmu',
      pvalueCutoff = pvalueSET,
      qvalueCutoff = qvalueSET
    )
    kegg.data <- setReadable(kegg_enrich, 'org.Mm.eg.db', 'ENTREZID')
  }
  
  kegg.data <- as.data.frame(kegg.data)
  if(select == 'pvalue'){
    kegg.data <- kegg.data[which(kegg.data$pvalue < cut_off_pvalue),]
    display_p <- 'pvalue'
  }else{
    kegg.data <- kegg.data[which(kegg.data$p.adjust < cut_off_pvalue),]
    display_p <- 'p.adjust'
  }
  
  rownames(kegg.data) <- 1:nrow(kegg.data)
  kegg.data$order <- factor(rev(as.integer(rownames(kegg.data))),labels = rev(kegg.data$Description))
  
  #气泡图
  p <- ggplot(data = kegg.data, aes(y = order, x = Count)) +
    geom_point(aes(size=Count,color=display_p)) + 
    scale_color_paletteer_c(
      palette   = "viridis::plasma",
      direction = -1,              # 让小 p 在高亮端
      name      = display_p
    ) +
    facet_grid(category ~., scales = "free_y",space = "free_y") +
    theme_bw() + # 黑白主题
    labs(x = "Gene Counts", y = "Pathways", title = "KEGG Enrichment") + # 设置x轴、y轴、标题标签
    theme(axis.title = element_text(size = 8), # 设置坐标轴标签字体大小
          axis.text = element_text(size = 8), # 设置坐标轴刻度字体大小
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"), # 设置标题字体大小、位置和加粗
          legend.title = element_text(size = 8), # 设置图例标题字体大小
          legend.text = element_text(size = 8), # 设置图例标签字体大小
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  results <- list(data = kegg.data, plot=p)
  return(results)

}

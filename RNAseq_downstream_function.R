#' @description   RNA-seq pipeline
#' @title  1）差异分析
#'   2）富集分析

#' @description
#' load & install required packages
#' @author shi jian
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
  species <- match.arg(species)
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
  mart <- useDataset(dataset = bio.use, useMart("ensembl"))
  gene2symbols <- getBM(attributes = c(id_type,'external_gene_name','entrezgene_id','description','gene_biotype'),
                        filters = id_type, values = gene_ids, mart = mart)
  
  return(gene2symbols)
}

#' @description 判断基因名在哪一列
#' @param df 输入数据
detect_gene_name_column <- function(df) {
  
  candidate_cols <- c("external_gene_name","gene_name","symbol","SYMBOL","Gene","gene")
  
  hit <- candidate_cols[candidate_cols %in% colnames(df)]
  
  if (length(hit) > 0) {
    message("✔ Use gene column: ", hit[1])
    return(hit[1])
  }
  
  ## fallback: rownames
  if (!is.null(rownames(df)) && all(rownames(df) != "")) {
    message("✔ Use rownames as gene names")
    return(NULL)  # NULL 表示用 rownames
  }
  
  stop("No gene name column detected, and rownames are empty.")
}

#' @description 检查数据类型
#' @param count_mat 表达矩阵 
#' @return 数据类型
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

#' @description 数据过滤
filter_fun <- function(counts,
                       min_count = 0,
                       min_samples = NULL,
                       mode = c("rowSums", "rowMeans")) {
  
  mode <- match.arg(mode)
  n_samples <- ncol(counts)
  
  if (is.null(min_samples)) {
    min_samples <- ceiling(n_samples / 2)
  }
  
  if (mode == "rowSums") {
    keep <- rowSums(counts >= min_count) >= min_samples
  } else {
    keep <- rowMeans(counts) >= min_count
  }
  
  return(keep)
}

#' @description 差异分析
#' @param counts 表达矩阵
#' @param coldata 样本信息
#' @param design 分组信息所在的表头名称
#' @param contrast 哪辆组比较
#' @param control_level 对照组
#' @param pvalue_type P值类型
run_DESeq2 <- function(
    counts,
    coldata,
    design,
    contrast = NULL,
    control_level = NULL,
    count.filter = 0,
    pvalue_type = c('padj','pvalue'),
    cut_off_pvalue = 0.05,
    cut_off_logFC = 1) {
  
  # import library
  nd <- c("DESeq2","dplyr","readxl",'tidyr')
  Plus.library(nd)
  
  pvalue_type <- match.arg(pvalue_type)
  
  # ---------- 1. input check ----------
  stopifnot(all(colnames(counts) %in% rownames(coldata)))
  coldata <- coldata[colnames(counts), , drop = FALSE] # 将coldata的行按照counts的列的顺序重新排序
  check_is_raw_counts(counts)
  
  # ---------- 2. relevel control ----------
  if (!is.null(control_level)) {
    var <- all.vars(design)[1]
    coldata[[var]] <- relevel(factor(coldata[[var]]), ref = control_level)
  }
  
  # ---------- 3. DESeq2 ----------
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = design
  )
  if (!is.null(count.filter)) {
    keep <- rowMeans(DESeq2::counts(dds)) >= count.filter
    dds <- dds[keep, ]
  }
  dds <- DESeq(dds)
  res <- if (is.null(contrast)) {
    results(dds)
  } else {
    results(dds, contrast = contrast)
  }
  
  # ---------- 4. DEG classification ----------
  deg <- as.data.frame(res) |>
    tibble::rownames_to_column("gene_id") |>
    mutate(
      change = case_when(
        .data[[pvalue_type]] < cut_off_pvalue & log2FoldChange >= cut_off_logFC  ~ "Up",
        .data[[pvalue_type]] < cut_off_pvalue & log2FoldChange <= -cut_off_logFC ~ "Down",
        TRUE ~ "Stable"
      )
    ) |>
    drop_na(.data[[pvalue_type]])
  
  return(list(dds = dds,result = res,DEG = deg))
}

#' @description 样本PCA分析
#' @param dds DESeq2对象
#' @param intgroup 分组信息（colData 中的列名）
#' @return 返回PCA数据和PCA图
plot_PCA <- function(dds,
                     intgroup = "condition"){
  nd <- c("ggplot2","ggrepel")
  Plus.library(nd)
  if (!intgroup %in% colnames(colData(dds))) {
    stop("intgroup not found in colData(dds)")
  }
  
  n_samples <- ncol(dds)
  
  if(n_samples < 30){
    scale.data <- rlog(dds, blind=TRUE)
  }else{
    scale.data <- vst(dds, blind=TRUE)
  }
  pcaData <- plotPCA(scale.data, intgroup= intgroup,returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca <- ggplot(pcaData, aes(PC1, PC2 , color = .data[[intgroup]],shape = .data[[intgroup]])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    ggtitle("Principal component plot")
  return(list(pcaData = pcaData,plot = pca))
}

#' @description 差异基因火山图
#' @param diff_data DESeq2结果
plot_volcanoplot <- function(diff_data,
                             cut_off_pvalue = 0.05,
                             cut_off_logFC = 1,
                             num = 5,
                             pvalue_type = c('padj','pvalue')){
  nd <- c("ggplot2","paletteer","tidyverse","dplyr")
  pvalue_type <- match.arg(pvalue_type)
  ## detect gene name
  gene_col <- detect_gene_name_column(diff_data)
  if (is.null(gene_col)) {
    diff_data$gene_label <- diff_data$gene_id
  } else {
    diff_data$gene_label <- diff_data$gene_id
  }
  # add gene labels in plot
  # 只考虑显著基因
  sig_data <- diff_data %>% 
    filter(.data[[pvalue_type]] < cut_off_pvalue)
  
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
  
  # genes to label
  to_label <- bind_rows(up, down) %>%
    mutate(label = gene_label)
  
  # merge label info (do NOT overwrite original data)
  diff_data1 <- diff_data %>%
    left_join(
      to_label %>% dplyr::select(gene_id, label),
      by = "gene_id"
    ) %>%
    mutate(
      label = replace_na(label, "")
    )
  
  p <- ggplot(
    diff_data1, aes(x = log2FoldChange, y = -log10(.data[[pvalue_type]]), colour=change)) +
    geom_point(alpha=0.8, size=0.5) +
    scale_color_manual(values=c("#104F8D", "#d2dae2","#BD4F48"))+
    geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
    labs(x="log2(Fold Change)",
         y=paste0("-log10 (",pvalue_type,")"))+
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
#' @description 差异基因热图
#' @param dds DESeq2对象，标准化数据时用到
#' @param coldata 样本信息
#' @param gene_list 基因列表
#' @param show_gene_name 是否展示基因名
plot_heatmap <- function(
    dds,
    coldata,
    gene_list,
    show_gene_name
) {
  nd <- c("DESeq2","pheatmap")
  Plus.library(nd)
  n_samples <- ncol(dds)
  if(n_samples < 30){
    scale.data <- rlog(dds, blind=TRUE)
  }else{
    scale.data <- vst(dds, blind=TRUE)
  }
  mat <- assay(scale.data)
  
  hit_genes <- intersect(rownames(mat), gene_list)
  if (length(hit_genes) == 0) {
    stop("No genes in gene_list matched rownames of dds.")
  }
  if (length(hit_genes) < length(gene_list)) {
    message("⚠ Only ", length(hit_genes), "/", length(gene_list),
            " genes matched and will be plotted.")
  }
  mat <- mat[hit_genes, , drop = FALSE]
  
  anno_col <- NULL
  if (!is.null(coldata)) {
    stopifnot(all(colnames(mat) %in% rownames(coldata)))
    anno_col <- coldata[colnames(mat), , drop = FALSE]
  }
  
  plot <- pheatmap(
    mat,
    scale = 'row',
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = anno_col,
    show_rownames = show_gene_name
    )
  
  return(plot)
}

GO_enrichment <- function(gene_list,
                          species = c("human", "mouse"),
                          pvalue_type = c('p.adjust','pvalue'),
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
    go <- enrichGO(gene = gene_list,
                   OrgDb = org.Hs.eg.db,     
                   keyType = "SYMBOL",     
                   ont = "ALL",              
                   pAdjustMethod = "BH",     
                   pvalueCutoff = pvalueSET, 
                   qvalueCutoff = qvalueSET             
    )
  }
  if(species == "mouse"){
    nd <- c("clusterProfiler","org.Mm.eg.db","ggplot2","paletteer")
    Plus.library(nd)
    go <- enrichGO(gene = gene_list,
                   OrgDb = org.Mm.eg.db,     
                   keyType = "SYMBOL",     
                   ont = "ALL",              
                   pAdjustMethod = "BH",     
                   pvalueCutoff = pvalueSET, 
                   qvalueCutoff = qvalueSET            
    )
  }
  
  go.data <- data.frame(go)
  if(pvalue_type == 'pvalue'){
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
    geom_point(aes(size = Count, color = .data[[pvalue_type]])) +
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
                          pvalue_type = c('p.adjust','pvalue'),
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
  if(pvalue_type == 'pvalue'){
    kegg.data <- kegg.data[which(kegg.data$pvalue < cut_off_pvalue),]
    display_p <- 'pvalue'
  }else{
    kegg.data <- kegg.data[which(kegg.data$p.adjust < cut_off_pvalue),]
    display_p <- 'p.adjust'
  }
  
  rownames(kegg.data) <- 1:nrow(kegg.data)
  kegg.data$order <- factor(rev(as.integer(rownames(kegg.data))),labels = rev(kegg.data$Description))
  
  kegg.data <- kegg.data[order(kegg.data$`pvalue`), ]
  # 从排序后的数据框中选择前10行 
  kegg.data.top <- head(kegg.data, n = showNum)
  #气泡图
  p <- ggplot(data = kegg.data.top, aes(y = order, x = Count)) +
    geom_point(aes(size=Count,color=.data[[pvalue_type]])) + 
    scale_color_paletteer_c(
      palette   = "viridis::plasma",
      direction = -1,              # 让小 p 在高亮端
      name      = display_p
    ) +
    #facet_grid(category ~., scales = "free_y",space = "free_y") +
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

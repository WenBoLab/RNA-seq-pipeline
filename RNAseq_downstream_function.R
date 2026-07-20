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

exchane_gene_name <- function(data,
                              species = c("human", "mouse"),
                              gene_info = NULL){
  if (!is.null(gene_info)) {
    
    colnames(gene_info) <- c("ensembl_gene_id","external_gene_name","gene_biotype")
    # 只保留输入的 gene_ids 对应的行
    gene_info <- gene_info[gene_info$ensembl_gene_id %in% rownames(data), , drop = FALSE]
    
    rownames(gene_info) <- gene_info$ensembl_gene_id
    exp <- merge(gene_info,data, by='row.names')
    is_duplicate <- duplicated(exp$external_gene_name)
    exp <- exp[!is_duplicate, ]
    rownames(exp) <- exp$external_gene_name
    exp <- exp[,-c(1:4)]
    return(exp)
  }
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
  id_type <- detect_gene_id_type(rownames(data))
  # gene2symbol
  mart <- useDataset(dataset = bio.use, useMart("ensembl"))
  gene2symbols <- getBM(attributes = c(id_type,'external_gene_name','gene_biotype'),
                        filters = id_type, values = rownames(data), mart = mart)
  rownames(gene2symbols) <- gene2symbols$ensembl_gene_id
  exp = merge(gene2symbols, data, by='row.names')
  exp$external_gene_name = make.unique(as.character(exp$external_gene_name), sep = ".")
  #is_duplicate <- duplicated(exp$external_gene_name)
  #exp <- exp[!is_duplicate, ]
  #rownames(exp) <- exp$external_gene_name
  exp <- exp[,-c(1:4)]
  return(exp)
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

#' @description Filter lowly expressed genes for bulk RNA-seq
#' @param data gene x sample raw count matrix
#' @param min_count minimum count to be considered expressed
#' @param min_samples minimum number of samples expressing the gene
#' @return filtered count matrix
filter_counts <- function(data,
                          min_count = 10,
                          min_samples = NULL) {
  
  stopifnot(is.matrix(data) || is.data.frame(data))
  
  n_samples <- ncol(data)
  
  if (is.null(min_samples)) {
    min_samples <- ceiling(n_samples / 2)
  }
  
  keep <- rowSums(data >= min_count) >= min_samples
  counts_filtered <- data[keep, , drop = FALSE]
  
  return(counts_filtered)
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
    count.filter = 10,
    sample.filter = 3,
    method = c('Wald','LRT'),
    pvalue_type = c('padj','pvalue'),
    cut_off_pvalue = 0.05,
    cut_off_logFC = 1) {
  
  # import library
  nd <- c("DESeq2","dplyr","readxl",'tidyr')
  Plus.library(nd)
  
  method <- match.arg(method)
  pvalue_type <- match.arg(pvalue_type)

  if (method == "LRT") {
    # LRT必须提供reduced_design，为空则终止报错
    if (is.null(reduced_design)) {
      stop("method='LRT' 时必须指定参数 reduced_design，Wald方法无需该参数！")
    }
  } else {
    # Wald模式，即使传了reduced_design也直接清空，避免干扰
    reduced_design <- NULL
  }
  
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
    keep <- rowSums(DESeq2::counts(dds) >= count.filter) >= sample.filter
    dds <- dds[keep, ]
  }
  if (method == 'LRT') {
    dds <- DESeq(dds,test = method, reduced = reduced_design)
  } else {
    dds <- DESeq(dds)
  }
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

#' @description 差异分析
#' @param counts 表达矩阵
#' @param coldata 样本信息
#' @param design 分组信息所在的表头名称
#' @param contrast 哪辆组比较
#' @param control_level 对照组
#' @param pvalue_type P值类型
run_limma <- function(
    counts,
    coldata,
    design,
    control_level = NULL,
    count.filter = 0,
    pvalue_type = c('adj.P.Val','P.Value'),
    cut_off_pvalue = 0.05,
    cut_off_logFC = 1) {
  
  # import library
  nd <- c("limma","dplyr","readxl",'tidyr')
  Plus.library(nd)
  
  pvalue_type <- match.arg(pvalue_type)
  
  # ---------- 1. input check ----------
  stopifnot(all(colnames(counts) %in% rownames(coldata)))
  coldata <- coldata[colnames(counts), , drop = FALSE] # 将coldata的行按照counts的列的顺序重新排序
  
  # ---------- 2. limma ----------
  group <- coldata[[design]]
  group <- gsub(" ", "_", group)
  group <- gsub("-", "_", group)
  group <- factor(group)
  
  if (!(control_level %in% levels(group))) {
    stop("control_level 不在分组中")
  }
  group <- relevel(group, ref = control_level)
  
  design_mat <- model.matrix(~ 0 + group)   # 去掉截距，便于显式对比
  colnames(design_mat) <- levels(factor(group))
  rownames(design_mat) <- colnames(counts)
  
  fit <- lmFit(counts, design_mat)
  
  groups <- setdiff(levels(group), control_level)
  res_list <- list()
  
  for (g in groups) {
    
    cont <- paste0(g, "-", control_level)
    
    contrast_vec <- makeContrasts(contrasts = cont, levels = design_mat)
    fit2 <- contrasts.fit(fit, contrast_vec)
    fit2 <- eBayes(fit2, trend = TRUE)
    
    deg <- topTable(fit2, number = Inf) |>
      tibble::rownames_to_column("gene_id") |>
      dplyr::mutate(
        comparison = cont,
        change = dplyr::case_when(
          .data[[pvalue_type]] < cut_off_pvalue & logFC >= cut_off_logFC  ~ "Up",
          .data[[pvalue_type]] < cut_off_pvalue & logFC <= -cut_off_logFC ~ "Down",
          TRUE ~ "Stable"
        )
      ) |>
      tidyr::drop_na(.data[[pvalue_type]])
    
    res_list[[cont]] <- deg
  }
  
  return(res_list)
}

#' @description 样本PCA分析
#' @param dds DESeq2对象
#' @param intgroup 分组信息（colData 中的列名）
#' @param sample_labels 是否显示样本标签，默认FALSE
#' @return 返回PCA数据和PCA图
plot_PCA <- function(data,
                     intgroup，
                     group_colors = FALSE,
                     sample_labels = FALSE){
  nd <- c("ggplot2","ggrepel")
  Plus.library(nd)
  # 判断输入类型
  is_deseq2 <- inherits(data, "DESeqDataSet")
  if (is_deseq2) {
    Plus.library("DESeq2")
    
    if (!intgroup %in% colnames(colData(data))) {
      stop("intgroup '", intgroup, "' not found in colData(dds)")
    }
    
    n_samples <- ncol(data)
    scale.data <- if (n_samples < 30) rlog(data, blind = TRUE) else vst(data, blind = TRUE)
    
    pcaData     <- plotPCA(scale.data, intgroup = intgroup, returnData = TRUE)
    percentVar  <- round(100 * attr(pcaData, "percentVar"))
    
    pca_df      <- pcaData
    group_col   <- intgroup
    pc1_label   <- paste0("PC1: ", percentVar[1], "% variance")
    pc2_label   <- paste0("PC2: ", percentVar[2], "% variance")
    pca_df$sample_id <- rownames(pca_df)
  } else{
    # ---- Array 分支 ----
    if (!is.matrix(data) && !is.data.frame(data)) {
      stop("Array模式下，data 请传入表达矩阵（行=基因，列=样本）")
    }
    expr_mat <- as.matrix(data)
    if (length(intgroup) == 1 && is.character(intgroup)) {
      stop("Array模式下，intgroup 请直接传入分组向量，例如：c('Control','Control','Treatment')")
    }
    if (length(intgroup) != ncol(expr_mat)) {
      stop("intgroup 长度 (", length(intgroup), ") 与样本数 (", ncol(expr_mat), ") 不一致")
    }
    pca_res    <- prcomp(t(expr_mat), scale. = TRUE, center = TRUE)
    var_exp    <- summary(pca_res)$importance[2, 1:2] * 100
    
    pca_df           <- as.data.frame(pca_res$x[, 1:2])
    pca_df$sample_id <- rownames(pca_df)
    pca_df$group     <- as.factor(intgroup)
    group_col        <- intgroup
    pc1_label        <- paste0("PC1: ", round(var_exp[1], 1), "% variance")
    pc2_label        <- paste0("PC2: ", round(var_exp[2], 1), "% variance")
  }
  
  # 绘图
  pca_plot <- ggplot(pca_df, aes(
    x     = PC1,
    y     = PC2,
    color = .data[[group_col]],
    #shape = .data[[group_col]],
    label = sample_id
  )) +
    stat_ellipse(aes(fill = .data[[group_col]]),
                 type = "norm", geom = "polygon",
                 alpha = 0.08, show.legend = FALSE) +
    geom_point(size = 2, alpha = 0.9) +
    { if (sample_labels)
      geom_text_repel(size = 3, color = "black",
                      box.padding = 0.4, max.overlaps = 20)
    } +
    xlab(pc1_label) +
    ylab(pc2_label) +
    ggtitle("Principal Component Analysis") +
    theme_classic(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.title    = element_blank(),
      legend.position = "right"
    )
  ## 如果用户指定颜色，则使用指定颜色
  if (!is.null(group_colors)) {
    pca_plot <- pca_plot +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors)
  }
  
  return(list(pcaData = pca_df,plot = pca_plot))
}

}

#' @description 差异基因火山图
#' @param diff_data DESeq2 或 limma 的差异分析结果数据框
#' @param cut_off_pvalue p值阈值
#' @param cut_off_logFC logFC阈值
#' @param num 标注的上/下调基因各前num个
#' @param pvalue_type 指定用哪列p值，NULL时自动检测
plot_volcanoplot <- function(diff_data,
                             cut_off_pvalue = 0.05,
                             cut_off_logFC = 1,
                             num = 5,
                             pvalue_type = c('padj','pvalue')){
  nd <- c("ggplot2","paletteer","tidyverse","dplyr")
  pvalue_type <- match.arg(pvalue_type)
  # 1. 统一字段名：logFC 和 pvalue
  ##logFC列
  logfc_col <- dplyr::case_when(
    "log2FoldChange" %in% colnames(diff_data) ~ "log2FoldChange",  # DESeq2
    "logFC"          %in% colnames(diff_data) ~ "logFC",           # limma
    TRUE ~ NA_character_
  )
  if (is.na(logfc_col)) stop("未找到 logFC 列，请确认输入为 DESeq2 或 limma 结果")
  
  ##p值列：优先用传入的，否则自动检测
  pvalue_candidates <- c("padj", "adj.P.Val", "pvalue", "P.Value")
  
  if (is.null(pvalue_type)) {
    pvalue_type <- pvalue_candidates[pvalue_candidates %in% colnames(diff_data)][1]
    if (is.na(pvalue_type)) stop("未找到 p值列，请手动指定 pvalue_type")
    message("自动检测 p值列: ", pvalue_type)
  } else {
    if (!pvalue_type %in% colnames(diff_data)) {
      stop("指定的 pvalue_type '", pvalue_type, "' 不在数据列中\n",
           "可用列: ", paste(colnames(diff_data), collapse = ", "))
    }
  }
  diff_data <- diff_data %>%
    dplyr::rename(
      logFC_unified  = !!logfc_col,
      pvalue_unified = !!pvalue_type
    )
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
    filter(pvalue_unified < cut_off_pvalue)
  
  # 上调（正 logFC）的前 num 个
  up <- sig_data %>%
    filter(logFC_unified > cut_off_logFC) %>%
    arrange(desc(logFC_unified)) %>%          # 从最大 logFC 开始
    slice_head(n = num) %>%
    mutate(direction = "Up")
  
  # 下调（负 logFC）的前 num 个
  down <- sig_data %>%
    filter(logFC_unified < -cut_off_logFC) %>%
    arrange(logFC_unified) %>%                # 从最小 logFC 开始
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
    diff_data1, aes(x = logFC_unified, y = -log10(pvalue_unified), colour=change)) +
    geom_point(alpha=0.8, size=0.5) +
    scale_color_manual(values=c("Down" = "#104F8D", "Stable" = "#d2dae2","Up" = "#BD4F48"))+
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
#' @param data DESeq2对象 或 归一化后的表达矩阵（行=基因，列=样本）
#' @param coldata 样本信息
#' @param gene_list 基因列表
#' @param anno_cols 指定coldata中用于注释的列名
#' @param show_gene_name 是否展示基因名
#' @param show_sample_name 是否展示样本名
#' @param ann_colors 分组颜色
#' @param cluster_cols 是否对列聚类
plot_heatmap <- function(
    data,
    coldata,
    gene_list,
    anno_cols = NULL,
    ann_colors = NULL,
    show_gene_name = TRUE,
    show_sample_name = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE
) {
  nd <- c("pheatmap")
  Plus.library(nd)
  
  # 1. 根据输入类型提取表达矩阵
  is_deseq2 <- inherits(data, "DESeqDataSet")
  if (is_deseq2) {
    # DESeq2 分支：rlog / vst 标准化
    Plus.library("DESeq2")
    n_samples  <- ncol(data)
    scale.data <- if (n_samples < 30) rlog(data, blind = FALSE) else vst(data, blind = FALSE)
    mat        <- assay(scale.data)
    
  } else {
    # Array / limma 分支：直接使用传入矩阵（应已归一化）
    if (!is.matrix(data) && !is.data.frame(data)) {
      stop("Array/limma 模式下，data 请传入表达矩阵（行=基因，列=样本）")
    }
    mat <- as.matrix(data)
    message("Array/limma 模式：直接使用传入矩阵，请确保已完成归一化")
  }
  
  # 2. 基因过滤
  hit_genes <- intersect(rownames(mat), gene_list)
  if (length(hit_genes) == 0) {
    stop("No genes in gene_list matched rownames of dds.")
  }
  if (length(hit_genes) < length(gene_list)) {
    message("⚠ Only ", length(hit_genes), "/", length(gene_list),
            " genes matched and will be plotted.")
  }
  mat <- mat[hit_genes, , drop = FALSE]
  # 3. 样本注释
  anno_col <- NULL
  if (!is.null(coldata)) {
    stopifnot(all(colnames(mat) %in% rownames(coldata)))
    anno_col <- coldata[colnames(mat), , drop = FALSE]
    ## 只保留指定列
    if (!is.null(anno_cols)) {
      stopifnot(all(anno_cols %in% colnames(anno_col)))
      anno_col <- anno_col[, anno_cols, drop = FALSE]
    }
  }
  
  # 4. 绘图
  # 行数较多时自动关闭基因名
  auto_show_gene <- if (missing(show_gene_name)) nrow(mat) <= 50 else show_gene_name
  plot <- pheatmap(
    mat,
    scale = 'row',
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_col = anno_col,
    show_rownames = auto_show_gene,
    show_colnames = show_sample_name,
    color         = colorRampPalette(c("#104F8D", "white", "#BD4F48"))(100),
    border_color   = NA,
    fontsize_row   = max(6, min(10, 200 / nrow(mat))),  # 字体大小随基因数自适应
    fontsize_col   = 10,
    annotation_colors = ann_colors
    )
  
  return(plot)
}


#' @description 富集分析结果排序
#' A short description...
#' 
sort_enrichment_result <- function(df,
                                   order_by = c("p.adjust","pvalue","GeneRatio","Count","RichFactor","NES"),
                                   decreasing = NULL){
  order_by <- match.arg(order_by)
  if(is.null(decreasing)){
    decreasing <- order_by %in%
      c("GeneRatio","Count","RichFactor")
  }
  
  ## GeneRatio需要转换
  if(order_by == "GeneRatio"){
    df$GeneRatio_num <- sapply(df$GeneRatio,function(x){
        tmp <- strsplit(as.character(x), "/")[[1]]
        as.numeric(tmp[1]) / as.numeric(tmp[2])
        
      })
    sort_column <- "GeneRatio_num"
  }else{
    sort_column <- order_by
  }
  if(!sort_column %in% colnames(df)){
    stop(paste0("排序参数 ",order_by," 不存在，可选:",paste(colnames(df),collapse=",")))
  }
  
  df <- df[order(df[[sort_column]],decreasing = decreasing),]
  return(df)
}

GO_enrichment <- function(gene_list,
                          species = c("human", "mouse"),
                          ont = c("ALL", "BP", "CC", "MF"),
                          pvalue_type = c('p.adjust','pvalue'),
                          cut_off_pvalue = 0.05,
                          showNum = 10,
                          order_by = c("p.adjust","pvalue","GeneRatio","Count","RichFactor"),
                          decreasing = NULL,
                          plot_type = c("dot", "bar")){
  
  # 参数匹配
  species <- match.arg(species)
  ont <- match.arg(ont)
  pvalue_type <- match.arg(pvalue_type)
  plot_type <- match.arg(plot_type)
  order_by <- match.arg(order_by)
  pvalueSET = 1
  qvalueSET = 1
  if(species == "human"){
    nd <- c("clusterProfiler","org.Hs.eg.db","ggplot2","paletteer")
    Plus.library(nd)
    go <- enrichGO(gene = gene_list,
                   OrgDb = org.Hs.eg.db,     
                   keyType = "SYMBOL",     
                   ont = ont,              
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
                   ont = ont,              
                   pAdjustMethod = "BH",     
                   pvalueCutoff = pvalueSET, 
                   qvalueCutoff = qvalueSET            
    )
  }
  
  go.data <- data.frame(go)
  if(nrow(go.data)==0){
    warning("No enriched GO terms found")
    return(NULL)
  }
  if(pvalue_type == 'pvalue'){
    go.data <- go.data[which(go.data$pvalue < cut_off_pvalue),]
    display_p <- 'pvalue'
  }else{
    go.data <- go.data[which(go.data$p.adjust < cut_off_pvalue),]
    display_p <- 'p.adjust'
  }
  if(nrow(go.data)==0){
    warning("No significant GO terms found")
    return(NULL)
  }
  
  if(ont != "ALL"){
    #go.data <- go.data[order(go.data[[pvalue_type]]), ]
    go.data <- sort_enrichment_result(go.data,order_by,decreasing)
    go.data.top <- head(go.data, n = showNum)
  } else{
    go.data.bp <- go.data[go.data$`ONTOLOGY` == "BP", ]
    go.data.cc <- go.data[go.data$`ONTOLOGY` == "CC", ]
    go.data.mf <- go.data[go.data$`ONTOLOGY` == "MF", ]
    
    #go.data.bp <- go.data.bp[order(go.data.bp[[pvalue_type]]), ]
    #go.data.cc <- go.data.cc[order(go.data.cc[[pvalue_type]]), ]
    #go.data.mf <- go.data.mf[order(go.data.mf[[pvalue_type]]), ]
    go.data.bp <- sort_enrichment_result(go.data.bp,order_by,decreasing)
    go.data.cc <- sort_enrichment_result(go.data.cc,order_by,decreasing)
    go.data.mf <- sort_enrichment_result(go.data.mf,order_by,decreasing)

    go.data.bp.top <- head(go.data.bp, n = showNum)
    go.data.cc.top <- head(go.data.cc, n = showNum)
    go.data.mf.top <- head(go.data.mf, n = showNum)
    
    go.data.top <- rbind(go.data.bp.top, go.data.cc.top, go.data.mf.top)
  }
  
  go.data.top$Description <- factor(go.data.top$Description,levels = rev(go.data.top$Description))
  
  gene_ratio_num <- sapply(strsplit(as.character(go.data.top$GeneRatio), "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]))
  
  go.data.top$GeneRatio_num <- gene_ratio_num
  
  if(plot_type == "dot"){
    p <- ggplot(go.data.top, aes(x = GeneRatio_num, y = Description)) +
      geom_point(aes(size = Count, color = .data[[order_by]])) +
      coord_cartesian(clip = "off") +
      # 颜色：用 viridis::plasma，并把颜色方向反过来（p越小颜色越亮）
      scale_color_paletteer_c(
        palette   = "viridis::plasma",
        direction = -1,              # 让小 p 在高亮端
        name      = order_by
      ) +
      scale_size_continuous(name = "Count",range = c(1, 5)) +
      labs(x = "GeneRatio",y = "Terms") +
      facet_grid(ONTOLOGY ~., scales = "free_y") +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 9)
      ) +
      scale_y_discrete(labels = scales::label_wrap(width = 60))
  } else{
    p <- ggplot(go.data.top, 
                aes(x = Count, 
                    #y = reorder(Description, Count))) +   # 兼顾P值和数量
                    y = Description)) +  # 只按照P值排序
      geom_col(aes(fill = .data[[order_by]]), 
               width = 0.75, color = "white", linewidth = 0.4) +
      scale_fill_gradient(low = "#ca0020", high = "#2166ac", 
                          name = order_by, trans = "log10") +
      labs(x = "Gene Count", y = "GO Terms", 
           title = paste("GO Enrichment", ont, "Bar Plot")) +
      theme_bw(base_size = 12) +
      theme(axis.text.y = element_text(size = 10.5),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "right")
  }
  
  results <- list(data = go.data, plot=p)
  return(results)
}
                           
KEGG_enrichment <- function(gene_list,
                          species = c("human", "mouse"),
                          pvalue_type = c('p.adjust','pvalue'),
                          cut_off_pvalue = 0.05,
                          showNum = 10,
                          order_by = c("p.adjust","pvalue","GeneRatio","Count","RichFactor"),
                          decreasing = NULL,
                          plot_type = c("bar", "dot")){
  
  # 参数匹配
  species <- match.arg(species)
  pvalue_type <- match.arg(pvalue_type)
  plot_type <- match.arg(plot_type)
  order_by <- match.arg(order_by)
  pvalueSET = 1
  qvalueSET = 1
  ## gene annotation
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
      gene = gene_anno$entrezgene_id,
      keyType = 'kegg',
      organism = 'mmu',
      pvalueCutoff = pvalueSET,
      qvalueCutoff = qvalueSET
    )
    kegg.data <- setReadable(kegg_enrich, 'org.Mm.eg.db', 'ENTREZID')
  }
  
  kegg.data <- as.data.frame(kegg.data)
  if(nrow(kegg.data)==0){
    warning("No enriched KEGG pathway")
    return(NULL)
  }
  
  if(pvalue_type == 'pvalue'){
    kegg.data <- kegg.data[which(kegg.data$pvalue < cut_off_pvalue),]
    display_p <- 'pvalue'
  }else{
    kegg.data <- kegg.data[which(kegg.data$p.adjust < cut_off_pvalue),]
    display_p <- 'p.adjust'
  }
  
  rownames(kegg.data) <- 1:nrow(kegg.data)
  #kegg.data <- kegg.data[order(kegg.data[[pvalue_type]]), ]
  kegg.data <- sort_enrichment_result(kegg.data,order_by,decreasing)

  # 从排序后的数据框中选择前10行 
  kegg.data.top <- head(kegg.data, n = showNum)
  #kegg.data.top$order <- factor(rev(as.integer(rownames(kegg.data.top))),labels = rev(kegg.data.top$Description))
  
  if(!"GeneRatio_num" %in%
     colnames(kegg.data.top)){
    kegg.data.top$GeneRatio_num <-
      sapply(
        strsplit(
          as.character(kegg.data.top$GeneRatio),
          "/"),
        function(x)
          as.numeric(x[1])/
          as.numeric(x[2])
      )
  }

  kegg.data.top$order <- factor(
    kegg.data.top$Description,
    levels = rev(kegg.data.top$Description)
  )
  
  if(plot_type == "dot") {
    # 气泡图
    p <- ggplot(data = kegg.data.top, aes(y = order, x = Count)) +
      geom_point(aes(size=Count,color=.data[[order_by]])) + 
      scale_color_paletteer_c(
        palette   = "viridis::plasma",
        direction = -1,              # 让小 p 在高亮端
        name      = order_by
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
  } else {
    # 柱状图
    p <- ggplot(data = kegg.data.top, 
                 aes(x = Count, y = order)) +     # 用 reorder 控制顺序
      geom_col(aes(fill = .data[[order_by]]), 
               width = 0.7, 
               color = "white",      # 柱子之间加白色间隔，更清晰
               linewidth = 0.3) +
      scale_fill_gradient(low = "#ca0020",high = "#2166ac",name = order_by,trans = "log10") +
      labs(x = "Gene Count", 
           y = "Pathways", 
           title = "KEGG Enrichment Analysis") +
      #facet_grid(category ~., scales = "free_y",space = "free_y") +
      theme_bw() +
      theme(
        axis.title = element_text(size = 11),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 9),
        plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "cm")
      ) +
      # 让 Y 轴标签不重叠
      theme(axis.text.y = element_text(lineheight = 0.8))
  }
  
  results <- list(data = kegg.data, plot=p)
  return(results)
}

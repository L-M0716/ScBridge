#!/usr/bin/env Rscript
library(optparse)
library(Seurat)
library(monocle3)
library(ggplot2)

# 解析命令行参数
option_list <- list(
  make_option(c("--input"), type = "character", help = "输入注释后的Seurat对象路径"),
  make_option(c("--output_trajectory"), type = "character", help = "输出拟时序对象路径"),
  make_option(c("--output_plot"), type = "character", help = "输出轨迹图PDF路径"),
  make_option(c("--root_markers"), type = "character", help = "根细胞标记基因（逗号分隔）"),
  make_option(c("--log"), type = "character", help = "日志文件路径")
)
opt <- parse_args(OptionParser(option_list = option_list))

# 日志函数
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = opt$log, append = TRUE)
}

# 主函数
main <- function() {
  tryCatch({
    # 加载数据
    log_msg(paste0("加载注释后的Seurat对象: ", opt$input))
    seurat_obj <- readRDS(opt$input)
    
    # 转换为CellDataSet对象
    log_msg("转换为Monocle3的CellDataSet对象...")
    cds <- as.cell_data_set(seurat_obj)
    
    # 估计大小因子
    log_msg("估计大小因子...")
    cds <- estimate_size_factors(cds)
    
    # 预处理
    log_msg("预处理数据...")
    cds <- preprocess_cds(cds, num_dim = 50)
    
    # 降维
    log_msg("执行降维...")
    cds <- reduce_dimension(cds)
    
    # 聚类
    log_msg("聚类细胞...")
    cds <- cluster_cells(cds)
    
    # 学习轨迹
    log_msg("学习拟时序轨迹...")
    cds <- learn_graph(cds)
    
    # 确定根细胞 - 修正：筛选表达根标记基因的细胞
    log_msg("确定根细胞...")
    root_markers <- unlist(strsplit(opt$root_markers, ","))
    
    # 修正：在细胞维度（列）上筛选
    root_cells <- colnames(seurat_obj)[colSums(GetAssayData(seurat_obj)[root_markers, ] > 0) > 0]
    
    # 如果没有找到根细胞，使用默认方法
    if (length(root_cells) == 0) {
      log_msg("警告：未找到表达根标记基因的细胞，使用默认根细胞选择方法")
      cds <- order_cells(cds)
    } else {
      cds <- order_cells(cds, root_cells = root_cells)
    }
    
    # 保存轨迹对象
    log_msg(paste0("保存拟时序对象: ", opt$output_trajectory))
    saveRDS(cds, file = opt$output_trajectory)
    
    # 生成轨迹图
    log_msg(paste0("生成轨迹图: ", opt$output_plot))
    pdf(opt$output_plot, width = 10, height = 8)
    
    p1 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
    p2 <- plot_cells(cds, color_cells_by = "manual_annotation", label_cell_groups = FALSE)
    
    print(p1)
    print(p2)
    dev.off()
    
    log_msg("拟时序分析完成")
    
  }, error = function(e) {
    log_msg(paste0("ERROR: ", e$message))
    stop(e$message)
  })
}

# 执行主函数
main()
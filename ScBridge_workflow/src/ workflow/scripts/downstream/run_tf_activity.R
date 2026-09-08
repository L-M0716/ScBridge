#!/usr/bin/env Rscript
library(optparse)
library(Seurat)
library(dorothea)
library(viper)
library(pheatmap)

# 解析命令行参数
option_list <- list(
  make_option(c("--input"), type = "character", help = "输入注释后的Seurat对象路径"),
  make_option(c("--output_table"), type = "character", help = "输出TF活性表路径"),
  make_option(c("--output_heatmap"), type = "character", help = "输出热图PDF路径"),
  make_option(c("--db"), type = "character", help = "转录因子数据库名称"),
  make_option(c("--db_path"), type = "character", help = "数据库路径"),
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
    
    # 加载转录因子数据库
    log_msg(paste0("加载转录因子数据库: ", opt$db_path))
    tf_db <- read.csv(opt$db_path)
    
    # 过滤高质量调控网络
    tf_db <- tf_db %>% 
      filter(confidence %in% c("A", "B"))  # 只保留高置信度
    
    # 计算TF活性 - 修正：使用viper函数
    log_msg("计算转录因子活性...")
    tf_activities <- viper(
      eset = GetAssayData(seurat_obj, slot = "scale.data"), 
      regulon = tf_db, 
      nes = TRUE
    )
    
    # 添加活性分数到Seurat对象
    seurat_obj[["tf_activity"]] <- CreateAssayObject(counts = tf_activities)
    
    # 按细胞类型计算平均TF活性
    log_msg("按细胞类型计算平均TF活性...")
    avg_tf_activity <- AverageExpression(seurat_obj, 
                                        assays = "tf_activity", 
                                        group.by = "manual_annotation")$tf_activity
    
    # 保存结果 - 修正：不保存行名
    log_msg(paste0("保存TF活性表: ", opt$output_table))
    write.csv(avg_tf_activity, file = opt$output_table, row.names = FALSE)  # 修正：不写行名
    
    # 生成热图
    log_msg(paste0("生成TF活性热图: ", opt$output_heatmap))
    pdf(opt$output_heatmap, width = 12, height = 16)
    
    # 选择高活性TF
    top_tfs <- names(sort(rowMeans(avg_tf_activity), decreasing = TRUE)[1:50])
    
    pheatmap(
      avg_tf_activity[top_tfs, ],
      main = "Transcription Factor Activity",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      fontsize_row = 8,
      color = colorRampPalette(c("blue", "white", "red"))(100)
    )
    dev.off()
    
    log_msg("转录因子活性分析完成")
    
  }, error = function(e) {
    log_msg(paste0("ERROR: ", e$message))
    stop(e$message)
  })
}

# 执行主函数
main()
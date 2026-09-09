#!/usr/bin/env Rscript
library(optparse)
library(Seurat)
library(CellChat)
library(patchwork)

# 解析命令行参数
option_list <- list(
  make_option(c("--input"), type = "character", help = "输入注释后的Seurat对象路径"),
  make_option(c("--output_obj"), type = "character", help = "输出CellChat对象路径"),
  make_option(c("--output_plot"), type = "character", help = "输出交互图PDF路径"),
  make_option(c("--db"), type = "character", help = "细胞通讯数据库名称"),
  make_option(c("--db_path"), type = "character", help = "数据库路径"),
  make_option(c("--interaction_categories"), type = "character", help = "交互类型（逗号分隔）"),
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
    
    # 创建CellChat对象
    log_msg(paste0("创建CellChat对象，使用数据库: ", opt$db))
    cellchat <- createCellChat(object = seurat_obj@assays$RNA@data, 
                              meta = seurat_obj@meta.data, 
                              group.by = "manual_annotation")
    
    # 加载细胞通讯数据库 - 修正：使用正确的物种名称
    log_msg(paste0("加载数据库: ", opt$db_path))
    species_name <- ifelse(opt$db == "mouse", "Mouse", "Human")
    CellChatDB <- loadCellChatDB(opt$db_path, species = species_name)  # 修正物种名称
    cellchat@DB <- CellChatDB
    
    # 预处理
    log_msg("预处理细胞通讯数据...")
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # 计算通讯概率
    log_msg("计算细胞通讯概率...")
    cellchat <- computeCommunProb(cellchat)
    
    # 过滤通讯 - 修正：使用配置的交互类型
    log_msg("过滤细胞通讯...")
    interaction_types <- unlist(strsplit(opt$interaction_categories, ","))
    cellchat <- filterCommunication(cellchat, 
                                   interaction.type = interaction_types,  # 使用配置的交互类型
                                   min.cells = 10)
    
    # 计算聚合通讯网络
    log_msg("计算聚合通讯网络...")
    cellchat <- computeCommunProbPathway(cellchat)
    
    # 聚合网络分析
    log_msg("聚合网络分析...")
    cellchat <- aggregateNet(cellchat)
    
    # 保存CellChat对象
    log_msg(paste0("保存CellChat对象: ", opt$output_obj))
    saveRDS(cellchat, file = opt$output_obj)
    
    # 生成交互图
    log_msg(paste0("生成交互图: ", opt$output_plot))
    pdf(opt$output_plot, width = 12, height = 10)
    
    # 1. 细胞通讯网络
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                     weight.scale = T, label.edge = F, 
                     title.name = "Number of interactions")
    
    # 2. 信号通路热图
    netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                      width = 15, height = 20)
    
    # 3. 特定通路可视化
    pathways <- cellchat@netP$pathways
    if (length(pathways) > 0) {
      for (pathway in pathways[1:min(3, length(pathways))]) {
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
      }
    }
    dev.off()
    
    log_msg("细胞通讯分析完成")
    
  }, error = function(e) {
    log_msg(paste0("ERROR: ", e$message))
    stop(e$message)
  })
}

# 执行主函数
main()
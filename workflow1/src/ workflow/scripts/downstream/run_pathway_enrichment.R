#!/usr/bin/env Rscript
library(optparse)
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠
library(org.Hs.eg.db)  # 人类
library(enrichplot)
library(ggplot2)
library(dplyr)

# 解析命令行参数
option_list <- list(
  make_option(c("--input"), type = "character", help = "输入差异分析结果CSV路径"),
  make_option(c("--output_table"), type = "character", help = "输出富集结果表路径"),
  make_option(c("--output_plot"), type = "character", help = "输出富集图PDF路径"),
  make_option(c("--organism"), type = "character", help = "生物体（mouse或human）"),
  make_option(c("--pval_cutoff"), type = "double", default = 0.05, help = "p值阈值"),
  make_option(c("--qval_cutoff"), type = "double", default = 0.2, help = "q值阈值"),
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
    log_msg(paste0("加载差异分析结果: ", opt$input))
    de_results <- read.csv(opt$input)
    
    # 获取显著差异基因 - 修正：使用配置的p值阈值
    sig_genes <- de_results %>% 
      filter(p_val_adj < opt$pval_cutoff) %>%  # 使用配置的阈值
      pull(gene)
    
    if (length(sig_genes) == 0) {
      stop("未发现显著差异基因，无法进行通路富集分析")
    }
    
    # 确定生物体数据库
    if (opt$organism == "mouse") {
      org_db <- org.Mm.eg.db
      kegg_org <- "mmu"
    } else if (opt$organism == "human") {
      org_db <- org.Hs.eg.db
      kegg_org <- "hsa"
    } else {
      stop(paste0("不支持的生物体: ", opt$organism))
    }
    
    # 转换基因ID
    log_msg("转换基因ID...")
    gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
    
    # 检查转换结果
    if (nrow(gene_ids) == 0) {
      stop("基因ID转换失败，未找到对应的ENTREZID")
    }
    
    # GO富集分析
    log_msg("执行GO富集分析...")
    go_enrich <- enrichGO(
      gene = gene_ids$ENTREZID,
      OrgDb = org_db,
      keyType = "ENTREZID",
      ont = "BP",
      pvalueCutoff = opt$pval_cutoff,
      qvalueCutoff = opt$qval_cutoff,
      readable = TRUE
    )
    
    # KEGG富集分析
    log_msg("执行KEGG富集分析...")
    kegg_enrich <- enrichKEGG(
      gene = gene_ids$ENTREZID,
      organism = kegg_org,
      keyType = "kegg",
      pvalueCutoff = opt$pval_cutoff,
      qvalueCutoff = opt$qval_cutoff
    )
    
    # 合并结果
    all_enrich <- rbind(
      as.data.frame(go_enrich) %>% mutate(source = "GO"),
      as.data.frame(kegg_enrich) %>% mutate(source = "KEGG")
    )
    
    # 保存结果
    log_msg(paste0("保存富集结果: ", opt$output_table))
    write.csv(all_enrich, file = opt$output_table, row.names = FALSE)
    
    # 生成富集图
    log_msg(paste0("生成富集图: ", opt$output_plot))
    pdf(opt$output_plot, width = 12, height = 10)
    
    if (nrow(all_enrich) > 0) {
      # 点图
      p1 <- dotplot(go_enrich, showCategory = 20) + 
        ggtitle("GO Biological Processes")
      
      # 网络图
      p2 <- cnetplot(go_enrich, categorySize = "pvalue", showCategory = 10)
      
      # KEGG通路图
      if (nrow(as.data.frame(kegg_enrich)) > 0) {
        p3 <- dotplot(kegg_enrich, showCategory = 15) + 
          ggtitle("KEGG Pathways")
        print(p1)
        print(p2)
        print(p3)
      } else {
        print(p1)
        print(p2)
      }
    } else {
      plot.new()
      text(0.5, 0.5, "No significant enrichment found", cex = 1.5)
    }
    dev.off()
    
    log_msg("通路富集分析完成")
    
  }, error = function(e) {
    log_msg(paste0("ERROR: ", e$message))
    stop(e$message)
  })
}

# 执行主函数
main()
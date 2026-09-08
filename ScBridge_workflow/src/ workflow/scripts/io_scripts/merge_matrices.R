        library(Matrix)
        library(Seurat)

        args <- commandArgs(trailingOnly = TRUE)
        output_dir <- args[1]
        matrix_files <- args[2:length(args)]

        # 读取所有矩阵
        matrices <- list()
        for (i in seq_along(matrix_files)) {
            dir_path <- dirname(matrix_files[i])
            mat <- Read10X(dir_path)
            # 为条形码添加后缀以区分来源
            colnames(mat) <- paste0(colnames(mat), "-", i)
            matrices[[i]] <- mat
        }

        # 合并矩阵
        merged_mat <- do.call(cbind, matrices)

        # 写出合并的矩阵
        writeMM(merged_mat, file.path(output_dir, "matrix.mtx"))
        write.table(colnames(merged_mat), file.path(output_dir, "barcodes.tsv"),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        write.table(data.frame(rownames(merged_mat), rownames(merged_mat), "Gene Expression"),
                    file.path(output_dir, "features.tsv"),
                    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "	")

        # 压缩文件
        system(paste("gzip -f", file.path(output_dir, "matrix.mtx")))
        system(paste("gzip -f", file.path(output_dir, "barcodes.tsv")))
        system(paste("gzip -f", file.path(output_dir, "features.tsv")))

        cat("合并完成
")
        cat("总细胞数:", ncol(merged_mat), "
")
        cat("总基因数:", nrow(merged_mat), "
")

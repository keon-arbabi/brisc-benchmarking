suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(BPCells)
})

options(future.globals.maxSize = Inf)
source(file.path(path.expand("~"), "sc-benchmarking", "utils_local.R"))

args = commandArgs(trailingOnly=TRUE)
DATA_NAME <- args[1]
DATA_PATH <- args[2]
OUTPUT_PATH_TIME <- args[3]
OUTPUT_PATH_DE <- args[4]
REP <- as.integer(args[length(args)])

bpcells_dir <- file.path(
  Sys.getenv("SCRATCH"), "bpcells", "de", paste0("data_", DATA_NAME))
unlink(bpcells_dir, recursive = TRUE)

system_info()
cat("--- Params ---\n")
cat("seurat de\n")
cat(sprintf("DATA_PATH=%s\n", DATA_PATH))

timers <- MemoryTimer(
  silent = FALSE, csv_path = OUTPUT_PATH_TIME,
  csv_columns = list(library = "seurat", test = "de",
                     dataset = DATA_NAME, rep = REP))

# BPCells native functions required:
# Loading from h5ad and writing to disk
timers$with_timer("Load data", {
  mat_disk <- open_matrix_anndata_hdf5(path = DATA_PATH)
  mat_disk <- convert_matrix_type(mat_disk, type = "uint32_t")
  write_matrix_dir(mat = mat_disk, dir = bpcells_dir)
  mat <- open_matrix_dir(dir = bpcells_dir)
  # Custom reader function for reading obs
  # without loading the entire h5ad
  obs_metadata <- read_h5ad_obs(DATA_PATH)
  data <- CreateSeuratObject(counts = mat, meta.data = obs_metadata)
})

rm(mat_disk, mat, obs_metadata); gc()

timers$with_timer("Quality control", {
  data[["percent.mt"]] <- PercentageFeatureSet(
    data, features = grep("^mt-", rownames(data),
    ignore.case = TRUE, value = TRUE))
  data[["malat1"]] <- FetchData(
    data, vars = grep("^malat1$", rownames(data),
    ignore.case = TRUE, value = TRUE))[, 1]
  data <- subset(
    data, subset = nFeature_RNA >= 100 & percent.mt <= 5 & malat1 > 0,
    slot = "counts")
})

if (DATA_NAME == "SEAAD") {
  data <- subset(data, subset = !is.na(cond))
  group_cols <- c("cond", "sample", "cell_type")
  ident_pairs <- list(test = "AD", ref = "Control")
} else if (DATA_NAME == "Parse") {
  data <- subset(data, subset = cytokine %in% c("IFN-gamma", "PBS"))
  group_cols <- c("cytokine", "sample", "cell_type")
  ident_pairs <- list(test = "IFN-gamma", ref = "PBS")
} else if (DATA_NAME == "PanSci") {
  data <- subset(data, subset = !is.na(cond))
  group_cols <- c("cond", "sample", "cell_type")
  ident_pairs <- list(test = "Aged", ref = "Young")
}

# Seurat lacks pseudobulk-level filtering

# BPCells native function required
# Seurat::AggregateExpression internally coerces the input
# count matrix through dgCMatrix which fails with
# x[i,j] too dense for [CR]sparseMatrix
# which exceeds R's 32-bit index limit
timers$with_timer("Pseudobulk", {
  cell_groups <- do.call(paste, c(data@meta.data[group_cols], sep = "_"))
  pb_mat <- GetAssayData(data, layer = "counts")
  pb_mat <- pseudobulk_matrix(pb_mat, cell_groups, method = "sum")
  first <- !duplicated(cell_groups)
  pb_meta <- data@meta.data[first, ]
  rownames(pb_meta) <- cell_groups[first]
  pb_meta <- pb_meta[colnames(pb_mat), ]
  data <- CreateSeuratObject(
    counts = as(pb_mat, "dgCMatrix"), meta.data = pb_meta)
  # Required for Seurat::FindMarkers internal filtering
  data <- NormalizeData(data, verbose = FALSE)
})

rm(cell_groups, pb_mat, pb_meta); gc()

timers$with_timer("Differential expression", {
  Idents(data) <- paste(
    data$cell_type, data@meta.data[[group_cols[1]]], sep = "_")
  de_list <- list()
  for (ct in unique(data$cell_type)) {
    de_list[[ct]] <- FindMarkers(
      data,
      ident.1 = paste(ct, ident_pairs$test, sep = "_"),
      ident.2 = paste(ct, ident_pairs$ref, sep = "_"),
      test.use = "DESeq2",
      verbose = FALSE)
  }
})

de_df <- do.call(rbind, lapply(names(de_list), function(ct) {
  df <- de_list[[ct]]
  data.frame(
    cell_type = ct,
    gene = rownames(df),
    logFC = df$avg_log2FC,
    p_value = df$p_val,
    p_value_adj = df$p_val_adj
  )
}))
write.csv(de_df, OUTPUT_PATH_DE, row.names = FALSE)

timers$shutdown()
cat("--- Completed successfully ---\n")

cat("\n--- Session Info ---\n")
print(sessionInfo())

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
})

options(future.globals.maxSize = Inf)
source(file.path(path.expand("~"), "brisc-benchmarking", "utils_local.R"))

args = commandArgs(trailingOnly=TRUE)
DATA_NAME <- args[1]
DATA_PATH <- args[2]
OUTPUT_PATH_TIME <- args[3]
OUTPUT_PATH_ACC <- args[4]
OUTPUT_PATH_LABELS <- args[5]
REP <- as.integer(args[length(args)])

bpcells_dir <- file.path(
  Sys.getenv("SCRATCH"), "bpcells", "transfer",
  paste0("data_", DATA_NAME, "_rep", REP))
unlink(bpcells_dir, recursive = TRUE)

system_info()

cat("--- Params ---\n")
cat("seurat transfer\n")
cat(sprintf("DATA_PATH=%s\n", DATA_PATH))

timers <- MemoryTimer(
  silent = FALSE, csv_path = OUTPUT_PATH_TIME,
  csv_columns = list(library = "seurat", test = "transfer",
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

timers$with_timer("Split data", {
  data_ref <- subset(data, subset = is_ref == 1)
  data_query <- subset(data, subset = is_ref == 0)
})

rm(data); gc()

timers$with_timer("Normalization", {
  data_ref <- NormalizeData(data_ref)
  data_query <- NormalizeData(data_query)
})

timers$with_timer("Feature selection", {
  data_ref <- FindVariableFeatures(data_ref)
})

timers$with_timer("PCA", {
  data_ref <- ScaleData(data_ref)
  data_ref <- RunPCA(data_ref)
})

cell_type_col <- if (DATA_NAME == "PanSci") "cell_type_fine" else "cell_type"

timers$with_timer("Transfer labels", {
    anchors <- FindTransferAnchors(
      reference = data_ref, query = data_query,
      reference.reduction = "pca")
    predictions <- TransferData(
      anchorset = anchors, refdata = data_ref@meta.data[[cell_type_col]])
    data_query <- AddMetaData(
      object = data_query, metadata = predictions)
})

accuracy_df <- transfer_accuracy(
  data_query@meta.data, cell_type_col, "predicted.id")
write.csv(accuracy_df, OUTPUT_PATH_ACC, row.names = FALSE)

labels_df <- data.frame(
  cell_id = colnames(data_query),
  cell_type = data_query@meta.data[[cell_type_col]],
  cell_type_transferred = data_query@meta.data[["predicted.id"]]
)
write.csv(labels_df, OUTPUT_PATH_LABELS, row.names = FALSE)

timers$shutdown()
cat("--- Completed successfully ---\n")

cat("\n--- Session Info ---\n")
print(sessionInfo())

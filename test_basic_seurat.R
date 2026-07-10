suppressPackageStartupMessages({
  library(BPCells)
  library(Seurat)
  library(presto)
})

options(future.globals.maxSize = Inf)
source(file.path(path.expand("~"), "sc-benchmarking", "utils_local.R"))

ARGS <- commandArgs(trailingOnly=TRUE)
DATA_NAME <- ARGS[1]
DATA_PATH <- ARGS[2]
OUTPUT_PATH_TIME <- ARGS[3]
OUTPUT_PATH_EMBEDDING <- ARGS[4]
OUTPUT_PATH_PCS <- ARGS[5]
OUTPUT_PATH_NEIGHBORS <- ARGS[6]
REP <- as.integer(ARGS[length(ARGS)])

system_info()
cat("--- Params ---\n")
cat("seurat basic\n")
cat(sprintf("DATA_PATH=%s\n", DATA_PATH))

timers <- MemoryTimer(
  silent = FALSE, csv_path = OUTPUT_PATH_TIME,
  csv_columns = list(library = "seurat", test = "basic",
                     dataset = DATA_NAME, rep = REP))

bpcells_dir <- file.path(
  Sys.getenv("SCRATCH"), "bpcells", "basic", paste0("data_", DATA_NAME))
unlink(bpcells_dir, recursive = TRUE)

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

timers$with_timer("Normalization", {
  data <- NormalizeData(data)
})

timers$with_timer("Feature selection", {
  data <- FindVariableFeatures(data)
})

timers$with_timer("PCA", {
  data <- ScaleData(data)
  data <- RunPCA(data)
})

# BPCells native functions required:
# Seurat::ComputeSNN fails with std::bad_alloc as the resulting sparse matrix
# exceeds R's 32-bit integer index limit.
timers$with_timer("Nearest neighbors", {
  pca <- Embeddings(data, "pca")
  knn <- knn_hnsw(pca, k = 21)
  snn <- knn_to_snn_graph(knn)
})

# BPCells native function required:
# Seurat::FindClusters (Louvain, 10 random starts, single-threaded)
# scales poorly as community count grows
timers$with_timer("Clustering", {
  for (resolution in c(0.25, 0.5, 1, 1.5, 2)) {
    cl <- cluster_graph_leiden(snn, resolution = resolution)
    data[[paste0("clusters_", format(resolution, nsmall = 1))]] <-
      as.character(cl)
  }
})

rm(snn); gc()

# RSpectra's ARPACK segfaults on poorly conditioned normalized Laplacians
# during the default spectral init. Patch uwot to use irlba instead since
# Seurat::RunUMAP does not expose uwot's init parameter
timers$with_timer("Embedding", {
  assignInNamespace("rspectra_is_installed", function() FALSE, ns = "uwot")
  data <- RunUMAP(data, dims = 1:50)
})

timers$with_timer("Find markers", {
  markers <- FindAllMarkers(data, group.by = "cell_type", only.pos = TRUE)
})

# Save PCs
pc_df <- as.data.frame(pca)
colnames(pc_df) <- paste0("PC_", seq_len(ncol(pca)))
arrow::write_parquet(pc_df, OUTPUT_PATH_PCS)

# Save neighbors
knn_idx <- remove_self_neighbors(knn$idx) - 1L
neighbors_df <- as.data.frame(knn_idx)
colnames(neighbors_df) <- paste0("neighbor_", seq_len(ncol(knn_idx)))
arrow::write_parquet(neighbors_df, OUTPUT_PATH_NEIGHBORS)

# Save embeddings
ct_col <- if (DATA_NAME == "PanSci") "cell_type_fine" else "cell_type"
embedding_df <- data.frame(
  cell_id = colnames(data),
  embed_1 = Embeddings(data, "umap")[, 1],
  embed_2 = Embeddings(data, "umap")[, 2],
  cell_type = data@meta.data[[ct_col]],
  cell_type_broad = data$cell_type_broad,
  cluster_res_0.25 = data$clusters_0.25,
  cluster_res_0.5 = data$clusters_0.5,
  cluster_res_1.0 = data$clusters_1.0,
  cluster_res_1.5 = data$clusters_1.5,
  cluster_res_2.0 = data$clusters_2.0
)
write.csv(embedding_df, OUTPUT_PATH_EMBEDDING, row.names = FALSE)

timers$shutdown()
cat("--- Completed successfully ---\n")

cat("\n--- Session Info ---\n")
print(sessionInfo())

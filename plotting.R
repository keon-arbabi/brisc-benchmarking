suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ragg)
})

setDTthreads(8)

work_dir <- file.path(path.expand("~"), "sc-benchmarking")

results <- list.files(
  file.path(work_dir, "output"), "_timer\\.csv$", full.names = TRUE) %>%
  discard(\(f) str_detect(basename(f), "brisc_.*_gpu_timer\\.csv$")) %>%
  map(\(f) fread(f) %>% mutate(.src = basename(f))) %>%
  rbindlist(fill = TRUE) %>% as_tibble()

grp <- list(
  order = c(
    "Seurat + BPCells",
    "Scanpy",
    "brisc (single-threaded)",
    "brisc (multi-threaded)"),
  colors = c(
    "brisc (multi-threaded)"  = "#e16032",
    "brisc (single-threaded)" = "#f38d70",
    "brisc (CPU)"             = "#e16032",
    "Scanpy"                  = "#269d69",
    "Seurat + BPCells"        = "#1c8ca8",
    "rapids-singlecell"       = "#7058be",
    "rapids-singlecell (GPU)" = "#7058be"))

grp_gpu <- c("rapids-singlecell (GPU)", "brisc (CPU)")

datasets <- c("SEAAD", "Parse", "PanSci")
ds_lab <- c(
  SEAAD = "1.2M cells (SEAAD)",
  Parse = "9.7M cells (Parse PBMC)",
  PanSci = "20.3M cells (PanSci)")
ds_unit <- c(SEAAD = "s", Parse = "h", PanSci = "m")
ds_div <- c(SEAAD = 1, Parse = 3600, PanSci = 60)

wf_lab <- c(basic = "Basic workflow", transfer = "Label transfer",
            de = "Differential expression")

op_order <- list(
  basic = c(
    "Load data", "Quality control", "Normalization",
    "Feature selection", "PCA", "Nearest neighbors",
    "Clustering", "Embedding", "Find markers"),
  transfer = c(
    "Load data", "Quality control", "Split data",
    "Normalization", "Feature selection", "PCA",
    "Transfer labels"),
  de = c(
    "Load data", "Quality control", "Pseudobulk",
    "Filter", "Differential expression"),
  commands = c(
    "Get expression by cell", "Get expression by gene",
    "Subset to one cell type",
    "Subset to highly variable genes",
    "Subsample to 10,000 cells",
    "Select categorical columns",
    "Split by cell type", "Concatenate cell types"))

op_lab <- setNames(unique(unlist(op_order)), unique(unlist(op_order)))
op_lab["Clustering"] <- "Cluster\n(5 resolutions)"
op_lab["Get expression by cell"] <- "Get expression\nby cell"
op_lab["Get expression by gene"] <- "Get expression\nby gene"
op_lab["Subset to one cell type"] <- "Subset to one\ncell type"
op_lab["Subset to highly variable genes"] <- "Subset to highly\nvariable genes"
op_lab["Subsample to 10,000 cells"] <- "Subsample to\n10,000 cells"
op_lab["Select categorical columns"] <- "Select categorical\ncolumns"
op_lab["Split by cell type"] <- "Split by\ncell type"
op_lab["Concatenate cell types"] <- "Concatenate\ncell types"

fmt_time <- function(s) {
  case_when(s >= 3600 ~ paste0(round(s / 3600, 1), "h"),
            s >= 60   ~ paste0(round(s / 60, 1), "m"),
            s >= 0.1  ~ paste0(round(s, 1), "s"),
            .default  = paste0(round(s * 1000), "ms"))
}

make_label <- function(time_s, ref_s) {
  lbl <- fmt_time(time_s)
  if_else(!is.na(ref_s) & ref_s > 0,
          paste0(lbl, " (", round(time_s / ref_s), "x)"), lbl)
}

base_theme <- theme_bw() + theme(
  text = element_text(color = "black"),
  axis.text = element_text(color = "black"),
  plot.title = element_text(size = 11, hjust = 0.5),
  panel.grid = element_blank(),
  plot.margin = margin(6, 6, 6, 6))

# --- Data preparation --------------------------------------------------------

prepared <- results %>%
  filter(library != "brisc" |
         !str_detect(operation, "^Embedding \\(") |
         operation == "Embedding (PaCMAP)") %>%
  mutate(
    operation = if_else(operation == "Embedding (PaCMAP)",
                        "Embedding", operation),
    workflow = recode(test, manipulation = "commands") %>%
      factor(levels = c(names(wf_lab), "commands")),
    dataset = factor(dataset, datasets),
    group = case_when(
      library == "brisc" & num_threads == -1 ~ "brisc (multi-threaded)",
      library == "brisc" ~ "brisc (single-threaded)",
      library == "scanpy" ~ "Scanpy",
      library == "rapids" ~ "rapids-singlecell",
      .default = "Seurat + BPCells") %>%
      factor(levels = c(grp$order, "rapids-singlecell", grp_gpu)),
    timed_out = FALSE)

# Impute timed-out steps for any (group, workflow, dataset) missing an op:
# first missing op gets duration = 24h - sum(completed), flagged timed_out.
TIMEOUT_S <- 24 * 3600

impute_timeouts <- function(d, wf) {
  d_wf <- d %>% filter(workflow == wf)
  if (nrow(d_wf) == 0) return(NULL)

  # Only consider ops that at least one group actually ran (avoids
  # imputing steps that no library performs).
  ops_expected <- intersect(op_order[[wf]], unique(d_wf$operation))

  existing <- d_wf %>% distinct(group, dataset, operation)
  all_combos <- d_wf %>%
    distinct(group, dataset) %>%
    crossing(operation = ops_expected)
  missing <- anti_join(all_combos, existing,
                       by = c("group", "dataset", "operation"))
  if (nrow(missing) == 0) return(NULL)

  first_missing <- missing %>%
    mutate(op_rank = match(operation, ops_expected)) %>%
    slice_min(op_rank, by = c(group, dataset), n = 1) %>%
    select(group, dataset, operation)

  sums <- d_wf %>%
    summarize(sum_dur = sum(duration, na.rm = TRUE),
              max_mem = max(memory, na.rm = TRUE),
              .by = c(group, dataset))

  first_missing %>%
    left_join(sums, by = c("group", "dataset")) %>%
    mutate(duration = pmax(TIMEOUT_S - sum_dur, 0),
           memory = max_mem,
           workflow = factor(wf, levels = levels(d$workflow)),
           test = wf, timed_out = TRUE,
           sum_dur = NULL, max_mem = NULL)
}

imputed <- map(names(wf_lab), \(wf) impute_timeouts(prepared, wf)) %>%
  list_rbind()
if (nrow(imputed) > 0) prepared <- bind_rows(prepared, imputed)

# --- Total runtime figure ----------------------------------------------------

prepared_main <- prepared %>%
  filter(group %in% grp$order) %>%
  mutate(group = factor(group, grp$order))

ref_total <- prepared_main %>%
  filter(workflow %in% names(wf_lab), group == "brisc (multi-threaded)") %>%
  summarize(ref = sum(duration), .by = c(workflow, dataset))

totals <- prepared_main %>%
  filter(workflow %in% names(wf_lab)) %>%
  summarize(total_s = sum(duration), .by = c(group, workflow, dataset)) %>%
  left_join(ref_total, by = c("workflow", "dataset")) %>%
  mutate(ref = if_else(group == "brisc (multi-threaded)", NA_real_, ref),
         dur_plot = total_s / ds_div[as.character(dataset)],
         label = make_label(total_s, ref))

peak_mem <- prepared_main %>%
  filter(workflow %in% names(wf_lab)) %>%
  summarize(memory = max(memory), .by = c(group, workflow, dataset))

plot_bar <- function(wf, ds) {
  totals %>%
    filter(workflow == wf, dataset == ds) %>%
    ggplot(aes(dur_plot, group, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label),
              hjust = -0.05, size = 2.5, na.rm = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_manual(values = grp$colors, guide = "none") +
    labs(title = if (wf == "basic") ds_lab[ds],
         x = if (wf == "de") paste0("Duration (", ds_unit[ds], ")"),
         y = if (ds == "SEAAD") wf_lab[wf]) +
    base_theme + theme(
      axis.title.y = element_text(size = 11, margin = margin(r = 10)),
      axis.text.y = if (ds != "SEAAD") element_blank()
                    else element_text(),
      axis.ticks.y = if (ds != "SEAAD") element_blank()
                     else element_line())
}

plot_mem <- function(wf, ds) {
  peak_mem %>%
    filter(workflow == wf, dataset == ds) %>%
    ggplot(aes(memory, group)) +
    geom_col(fill = "#bdbdbd", width = 0.7) +
    geom_text(aes(x = 0, label = paste0(round(memory, 1), " GiB")),
              hjust = -0.05, size = 2.5) +
    scale_x_continuous(breaks = scales::breaks_pretty(3),
                       expand = expansion(mult = c(0, 0.5))) +
    scale_y_discrete(drop = FALSE) +
    labs(title = if (wf == "basic") "Peak memory",
         x = if (wf == "de") "GiB", y = NULL) +
    base_theme + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

fig_totals <- map(names(wf_lab), \(wf) {
  map(datasets, \(ds) list(plot_bar(wf, ds), plot_mem(wf, ds))) %>%
    list_flatten()
}) %>% list_flatten() %>%
  wrap_plots(ncol = 6, byrow = TRUE, widths = rep(c(4, 1), 3))

ggsave(file.path(work_dir, "figures", "fig_total_runtime.png"),
       fig_totals, width = 15, height = 6, device = agg_png)

# --- Per-step figures --------------------------------------------------------

compute_steps <- function(d, ref_group) {
  step_ref <- d %>%
    filter(group == ref_group) %>%
    summarize(ref = mean(duration),
              .by = c(operation, workflow, dataset))
  d %>%
    summarize(duration = mean(duration), memory = max(memory),
              timed_out = any(timed_out),
              .by = c(group, operation, workflow, dataset)) %>%
    left_join(step_ref, by = c("operation", "workflow", "dataset")) %>%
    mutate(ref = if_else(group == ref_group, NA_real_, ref),
           dur_plot = duration / ds_div[as.character(dataset)],
           label = make_label(duration, ref),
           label = if_else(timed_out,
                           paste0(label, ", timed out"), label)) %>%
    mutate(dur_plot = pmax(dur_plot,
                           max(dur_plot, na.rm = TRUE) * 0.02),
           .by = c(operation, workflow, dataset))
}

gpu_rename <- c("rapids-singlecell"      = "rapids-singlecell (GPU)",
                "brisc (multi-threaded)" = "brisc (CPU)")
prepared_gpu <- prepared %>%
  filter(group %in% names(gpu_rename)) %>%
  mutate(group = factor(gpu_rename[as.character(group)], levels = grp_gpu))

steps_main <- compute_steps(prepared_main, "brisc (multi-threaded)")
steps_gpu  <- compute_steps(prepared_gpu,  "brisc (CPU)")

build_step_fig <- function(wf, steps_df, group_order, ref_label,
                           legend_title, include_memory = TRUE) {
  ops <- intersect(op_order[[wf]],
                   unique(steps_df$operation[steps_df$workflow == wf]))
  if (length(ops) == 0) return(NULL)

  d <- steps_df %>%
    filter(workflow == wf, operation %in% ops) %>%
    mutate(operation = factor(operation, ops, labels = op_lab[ops]))

  x_title <- paste0("Per-step runtime\n(\u00d7 slower than ", ref_label, ")")

  make_title <- function(ds) {
    ggplot() + theme_void() +
      labs(title = ds_lab[ds]) +
      theme(plot.title = element_text(hjust = 0.5, size = 11,
                                      margin = margin(0, 0, 0, 0)),
            plot.margin = margin(0, 0, 0, 0))
  }

  make_t <- function(j) {
    ds <- datasets[j]; is_first <- j == 1; is_last <- j == length(datasets)
    d_ds <- d %>% filter(dataset == ds)
    ggplot(d_ds, aes(dur_plot, group, fill = group)) +
      geom_col(width = 0.8) +
      geom_segment(
        data = ~ filter(.x, timed_out),
        aes(x = dur_plot * 0.98, xend = dur_plot * 1.08,
            y = group, yend = group, color = group),
        arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
        linewidth = 0.8,
        inherit.aes = FALSE, na.rm = TRUE, show.legend = FALSE) +
      scale_color_manual(values = grp$colors, guide = "none") +
      geom_text(aes(label = label),
                hjust = -0.15, size = 2, na.rm = TRUE) +
      facet_wrap(~ operation, ncol = 1, scales = "free_x",
                 strip.position = "left") +
      scale_x_continuous(expand = expansion(mult = c(0, 1.1))) +
      scale_y_discrete(drop = FALSE) +
      scale_fill_manual(
        values = grp$colors, limits = group_order, name = legend_title,
        guide = if (is_first) guide_legend(reverse = TRUE)
                else "none") +
      labs(x = x_title, y = NULL) +
      base_theme +
      theme(
        strip.text.y.left = if (is_first)
          element_text(angle = 0, hjust = 1, size = 8)
          else element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 8,
                                    margin = margin(t = 5)),
        panel.spacing.y = unit(-0.5, "pt"),
        plot.margin = margin(
          0,
          if (include_memory) 0 else if (is_last) 4 else 8,
          2,
          if (is_first) 4 else 8))
  }

  make_m <- function(j) {
    ds <- datasets[j]; is_last <- j == length(datasets)
    d_ds <- d %>% filter(dataset == ds)
    mem_max <- max(d_ds$memory, na.rm = TRUE)
    mem_nudge <- mem_max * 0.08
    ggplot(d_ds, aes(memory, group)) +
      geom_col(fill = "#d9d9d9", width = 0.8) +
      geom_text(aes(x = 0, label = paste0(round(memory), " GiB")),
                hjust = 0, nudge_x = mem_nudge, size = 2) +
      facet_wrap(~ operation, ncol = 1, scales = "free_x",
                 strip.position = "right") +
      scale_x_continuous(limits = c(0, mem_max * 1.5),
                         expand = expansion(mult = c(0, 0))) +
      scale_y_discrete(drop = FALSE) +
      labs(x = "Peak\nmemory", y = NULL) +
      base_theme +
      theme(strip.text = element_blank(),
            strip.background = element_blank(),
            strip.placement = "outside",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.length.y = unit(0, "pt"),
            axis.title.x = element_text(size = 8,
                                        margin = margin(t = 5)),
            panel.spacing.y = unit(-0.5, "pt"),
            plot.margin = margin(0, if (is_last) 4 else 0, 2, -1))
  }

  panels <- list()
  for (j in seq_along(datasets)) {
    panels[[paste0("title", j)]] <- make_title(datasets[j])
    panels[[paste0("t", j)]] <- make_t(j)
    if (include_memory) panels[[paste0("m", j)]] <- make_m(j)
  }

  if (include_memory) {
    design <- "AABBCC\nDEFGHI"
    layout <- list(
      A = panels$title1, B = panels$title2, C = panels$title3,
      D = panels$t1, E = panels$m1,
      F = panels$t2, G = panels$m2,
      H = panels$t3, I = panels$m3)
    widths <- c(4, 1.3, 4, 1.3, 4, 1.3)
  } else {
    design <- "ABC\nDEF"
    layout <- list(
      A = panels$title1, B = panels$title2, C = panels$title3,
      D = panels$t1, E = panels$t2, F = panels$t3)
    widths <- c(4, 4, 4)
  }

  wrap_plots(layout, design = design,
             heights = c(0.15, 10), widths = widths) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.title = element_text(size = 9),
          legend.text  = element_text(size = 9))
}

walk(c(names(wf_lab), "commands"), \(wf) {
  n <- length(intersect(
    op_order[[wf]],
    unique(steps_main$operation[steps_main$workflow == wf])))
  ggsave(file.path(work_dir, "figures",
                   paste0("fig_", wf, "_steps.png")),
         build_step_fig(wf, steps_main, grp$order, "brisc MT",
                        legend_title = "192 CPUs, 755 GB RAM"),
         width = 11, height = max(3, n * 0.85),
         device = agg_png)
})

n_gpu <- length(intersect(
  op_order$basic,
  unique(steps_gpu$operation[steps_gpu$workflow == "basic"])))
ggsave(file.path(work_dir, "figures", "fig_basic_gpu_steps.png"),
       build_step_fig("basic", steps_gpu, grp_gpu, "brisc CPU",
                      legend_title = "96 CPUs,\n4\u00d7 H100 GPU,\n752 GB RAM",
                      include_memory = FALSE),
       width = 10.5, height = max(2, n_gpu * 0.42),
       device = agg_png)

# --- Embedding figures --------------------------------------------------------
# Top row: brisc (PaCMAP, LocalMAP, UMAP, UMAP hogwild)
# Bottom row: Scanpy UMAP, Seurat + BPCells UMAP, rapids-singlecell UMAP

embed_paths <- c(
  brisc  = "basic_brisc_%s_-1_embedding.csv",
  scanpy = "basic_scanpy_%s_embedding.csv",
  seurat = "basic_seurat_%s_embedding.csv",
  rapids = "basic_rapids_%s_gpu_embedding.csv")

lib_display <- c(
  brisc = "brisc", scanpy = "Scanpy",
  seurat = "Seurat + BPCells", rapids = "rapids-singlecell")

panel_spec <- tribble(
  ~lib,     ~method,          ~col1,            ~col2,            ~op,
  "brisc",  "PaCMAP",         "pacmap_1",       "pacmap_2",       "Embedding (PaCMAP)",
  "brisc",  "LocalMAP",       "localmap_1",     "localmap_2",     "Embedding (LocalMAP)",
  "brisc",  "UMAP",           "umap_1",         "umap_2",         "Embedding (UMAP)",
  "brisc",  "UMAP (hogwild)", "umap_hogwild_1", "umap_hogwild_2", "Embedding (UMAP hogwild)",
  "scanpy", "UMAP",           "embed_1",        "embed_2",        "Embedding",
  "seurat", "UMAP",           "embed_1",        "embed_2",        "Embedding",
  "rapids", "UMAP",           "embed_1",        "embed_2",        "Embedding")

get_embed_time <- function(lib_val, op_val, ds_val) {
  t <- results %>%
    filter(test == "basic", dataset == ds_val,
           library == lib_val, operation == op_val,
           !(library == "brisc" & num_threads != -1),
           !(library == "brisc" & str_detect(.src, "_gpu_timer"))) %>%
    pull(duration)
  if (length(t) > 0) t[1] else NA_real_
}

for (ds in datasets) {
  files <- file.path(work_dir, "output", sprintf(embed_paths, ds))
  names(files) <- names(embed_paths)
  files <- files[file.exists(files)]
  if (length(files) == 0) next

  raw_data <- imap(files, \(f, lib) {
    fread(f) %>% as_tibble() %>%
      slice_sample(n = min(nrow(.), 50000))
  })

  spec <- panel_spec %>% filter(lib %in% names(raw_data))

  embed_data <- spec %>%
    pmap(function(lib, method, col1, col2, op) {
      df <- raw_data[[lib]]
      if (!all(c(col1, col2) %in% names(df))) return(NULL)
      df %>%
        transmute(embed_1 = .data[[col1]], embed_2 = .data[[col2]],
                  cluster_res_1.0, cell_type_broad,
                  lib = !!lib, method = !!method,
                  panel_key = paste0(!!lib, "_", !!method))
    }) %>% bind_rows()
  if (nrow(embed_data) == 0) next

  embed_data <- embed_data %>%
    mutate(cluster_id = paste0(panel_key, "_", cluster_res_1.0))

  panel_info <- spec %>%
    mutate(panel_key = paste0(lib, "_", method),
           duration = map2_dbl(lib, op, get_embed_time, ds),
           label = paste0(
             lib_display[lib], " (", method,
             if_else(is.na(duration), "",
                     paste0(", ", fmt_time(duration))),
             ")"))

  embed_data <- embed_data %>%
    mutate(panel_label = factor(panel_key,
                                 levels = panel_info$panel_key,
                                 labels = panel_info$label))

  broad_types <- sort(unique(embed_data$cell_type_broad))
  hues <- setNames(
    head(seq(0, 360, length.out = length(broad_types) + 1), -1),
    broad_types)

  pal_df <- embed_data %>%
    count(panel_key, cluster_res_1.0, cell_type_broad, sort = TRUE) %>%
    distinct(panel_key, cluster_res_1.0, .keep_all = TRUE) %>%
    arrange(as.numeric(cluster_res_1.0)) %>%
    mutate(rank = row_number(), n_cl = n(),
           .by = c(panel_key, cell_type_broad)) %>%
    mutate(color = hcl(
      h = hues[cell_type_broad], c = 55,
      l = 85 - (85 - 45) * (rank - 1) / pmax(n_cl - 1, 1)))
  pal <- setNames(pal_df$color,
                  paste0(pal_df$panel_key, "_", pal_df$cluster_res_1.0))

  broad_pal <- setNames(
    hcl(h = hues[broad_types], c = 55, l = 60), broad_types)

  base_plot <- function(d) {
    ggplot(d, aes(embed_1, embed_2)) +
      geom_point(aes(color = cluster_id),
                 size = 0.375, stroke = 0, alpha = 1) +
      scale_color_manual(values = pal, guide = "none") +
      geom_point(aes(fill = cell_type_broad),
                 shape = 22, size = 0, stroke = 0) +
      scale_fill_manual(values = broad_pal, name = NULL,
        guide = guide_legend(override.aes = list(size = 4))) +
      facet_wrap(~ panel_label, nrow = 1, scales = "free") +
      labs(x = NULL, y = NULL) +
      base_theme +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            strip.text = element_text(size = 10),
            strip.background = element_blank(),
            legend.position = "right")
  }

  top_keys <- panel_info$panel_key[panel_info$lib == "brisc"]
  bot_keys <- panel_info$panel_key[panel_info$lib != "brisc"]

  p_top <- base_plot(embed_data %>% filter(panel_key %in% top_keys)) +
    theme(legend.position = "none")
  p_bot <- base_plot(embed_data %>% filter(panel_key %in% bot_keys))

  n_top <- length(top_keys)
  n_bot <- length(bot_keys)
  pad <- (n_top - n_bot) / 2
  p_bot_centered <- if (n_bot < n_top) {
    plot_spacer() + p_bot + plot_spacer() +
      plot_layout(widths = c(pad, n_bot, pad))
  } else p_bot

  p <- (p_top / p_bot_centered) +
    plot_annotation(title = ds_lab[ds],
                    theme = theme(plot.title = element_text(
                      size = 11, hjust = 0.5)))

  ggsave(file.path(work_dir, "figures",
                   paste0("fig_embeddings_", ds, ".png")),
         p, width = 12, height = 6, dpi = 300, device = agg_png)
}

# --- Accuracy figures --------------------------------------------------------
# Label-transfer accuracy (total and per cell type) and kNN recall.

acc_files <- list.files(
  file.path(work_dir, "output"),
  "^transfer_.*_accuracy\\.csv$", full.names = TRUE)

accuracy_data <- map(acc_files, \(f) {
  fread(f) %>% as_tibble() %>% mutate(.src = basename(f))
}) %>% bind_rows() %>%
  mutate(
    library = str_extract(.src, "(?<=^transfer_)[^_]+"),
    dataset = str_extract(.src, "SEAAD|Parse|PanSci"),
    nt = suppressWarnings(as.integer(
      str_extract(.src, "(?<=_)(-?\\d+)(?=_accuracy)"))),
    dataset = factor(dataset, datasets),
    group = case_when(
      library == "brisc" & nt == -1 ~ "brisc (multi-threaded)",
      library == "brisc" & nt == 1 ~ "brisc (single-threaded)",
      library == "scanpy" ~ "Scanpy",
      library == "seurat" ~ "Seurat + BPCells",
      .default = NA_character_)) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group, grp$order))

for (ds in datasets) {
  d_ds <- accuracy_data %>% filter(dataset == ds)
  if (nrow(d_ds) == 0) next

  # Overall accuracy per group (aggregated across cell types)
  total_ds <- d_ds %>%
    summarize(percent_correct = sum(n_correct) / sum(n_total) * 100,
              n_correct = sum(n_correct), n_total = sum(n_total),
              .by = c(group, dataset)) %>%
    mutate(cell_type = "Overall")

  d_with_total <- bind_rows(total_ds, d_ds)

  # Sort cell types by mean accuracy descending (highest first)
  ct_sorted <- d_ds %>%
    summarize(mean_pct = mean(percent_correct, na.rm = TRUE),
              .by = cell_type) %>%
    arrange(desc(mean_pct)) %>%
    pull(cell_type)

  make_col <- function(types, show_leg) {
    # Factor levels ascending => highest at top of y-axis
    d <- d_with_total %>%
      filter(cell_type %in% types) %>%
      mutate(cell_type = factor(cell_type, rev(types)))
    ggplot(d, aes(percent_correct, cell_type, fill = group)) +
      geom_col(position = position_dodge(width = 0.8,
                                         preserve = "single"),
               width = 0.7) +
      scale_x_continuous(limits = c(0, 100),
                         expand = expansion(mult = c(0, 0.02))) +
      scale_fill_manual(values = grp$colors, limits = grp$order,
                        name = NULL,
                        guide = if (show_leg) guide_legend(reverse = TRUE)
                                else "none") +
      labs(x = "Label-transfer accuracy (%)", y = NULL) +
      base_theme +
      theme(axis.text.y = element_text(size = 8),
            legend.position = "right")
  }

  if (ds == "PanSci") {
    half <- ceiling(length(ct_sorted) / 2)
    col1_types <- c("Overall", ct_sorted[seq_len(half)])        # higher
    col2_types <- ct_sorted[(half + 1):length(ct_sorted)]       # lower
    p <- (make_col(col1_types, FALSE) | make_col(col2_types, TRUE)) +
      plot_annotation(title = ds_lab[ds],
                      theme = theme(plot.title = element_text(
                        size = 11, hjust = 0.5)))
    h <- max(3, (half + 1) * 0.32)
    w <- 14
  } else {
    all_types <- c("Overall", ct_sorted)
    p <- make_col(all_types, TRUE) +
      labs(title = ds_lab[ds]) +
      theme(plot.title = element_text(size = 11, hjust = 0.5))
    h <- max(3, length(all_types) * 0.32)
    w <- 9
  }

  ggsave(file.path(work_dir, "figures",
                   paste0("fig_transfer_celltype_", ds, ".png")),
         p, width = w, height = h, device = agg_png)
}

# --- Clustering quality figures ---------------------------------------------
# kNN recall and Leiden cluster ARI vs ground-truth labels.

knn_colors <- c(
  "Seurat + BPCells"  = "#1c8ca8",
  "Scanpy"            = "#269d69",
  "brisc"             = "#e16032",
  "rapids-singlecell" = "#7058be")
knn_order <- names(knn_colors)
method_to_group <- c(seurat = "Seurat + BPCells", scanpy = "Scanpy",
                     brisc = "brisc", rapids = "rapids-singlecell")

knn_path <- file.path(work_dir, "output", "knn_recall_summary.csv")
if (file.exists(knn_path)) {
  knn_data <- fread(knn_path) %>%
    as_tibble() %>%
    mutate(dataset = factor(dataset, datasets),
           group = factor(method_to_group[method], levels = knn_order))

  plot_knn <- ggplot(knn_data, aes(mean, group, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", mean)),
              hjust = -0.15, size = 2.5, na.rm = TRUE) +
    facet_wrap(~ dataset, nrow = 1, labeller = as_labeller(ds_lab)) +
    scale_x_continuous(limits = c(0, 1.15),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_manual(values = knn_colors, guide = "none") +
    labs(x = "Mean kNN recall (vs exact kNN)", y = NULL) +
    base_theme +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 10))

  ggsave(file.path(work_dir, "figures", "fig_knn_recall.png"),
         plot_knn, width = 11, height = 2.5, device = agg_png)
}

ari_path <- file.path(work_dir, "output", "cluster_ari_summary.csv")
if (file.exists(ari_path)) {
  ari_data <- fread(ari_path) %>%
    as_tibble() %>%
    mutate(
      dataset = factor(dataset, datasets),
      group = factor(method_to_group[method], levels = knn_order),
      label = factor(label, c("cell_type_broad", "cell_type"),
                     labels = c("Broad cell type", "Fine cell type")))

  plot_ari <- ggplot(ari_data,
                     aes(resolution, balanced_ari, color = group)) +
    geom_errorbar(aes(ymin = balanced_ari - balanced_ari_std,
                      ymax = balanced_ari + balanced_ari_std),
                  width = 0.04, linewidth = 0.4) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.6) +
    facet_grid(label ~ dataset, labeller = labeller(dataset = ds_lab)) +
    scale_color_manual(values = knn_colors, limits = knn_order, name = NULL) +
    scale_x_continuous(breaks = c(0.25, 0.5, 1, 1.5, 2)) +
    scale_y_continuous(limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Leiden resolution",
         y = "Balanced ARI vs ground-truth cell types") +
    base_theme +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10),
          legend.position = "right",
          panel.grid.major.y = element_line(color = "grey92"))

  ggsave(file.path(work_dir, "figures", "fig_cluster_ari.png"),
         plot_ari, width = 11, height = 5, device = agg_png)
}

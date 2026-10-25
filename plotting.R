suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ragg)
  library(svglite)
  library(scattermore)
})

setDTthreads(8)

work_dir <- file.path(path.expand("~"), "brisc-benchmarking")
fig_dir <- file.path(work_dir, "figures")
font <- "Arial"

datasets <- c("SEAAD", "Parse", "PanSci")
ds_lab <- c(SEAAD = "1.4M cells (SEA-AD)",
            Parse = "9.7M cells (Parse)",
            PanSci = "20.3M cells (PanSci)")
ds_unit <- c(SEAAD = "s", Parse = "h", PanSci = "m")
ds_div <- c(SEAAD = 1, Parse = 3600, PanSci = 60)

wf_lab <- c(basic = "Standard workflow", transfer = "Label transfer",
            de = "Differential expression")

grp_order <- c("Seurat + BPCells", "Scanpy", "brisc (single-threaded)",
               "brisc (multi-threaded)")
grp_gpu <- c("rapids-singlecell (GPU)", "brisc (CPU)")
grp_colors <- c("brisc (multi-threaded)"  = "#e16032",
                "brisc (single-threaded)" = "#f38d70",
                "brisc (CPU)"             = "#e16032",
                "Scanpy"                  = "#269d69",
                "Seurat + BPCells"        = "#1c8ca8",
                "rapids-singlecell"       = "#7058be",
                "rapids-singlecell (GPU)" = "#7058be")

lib_colors <- c("brisc"             = "#e16032",
                "Scanpy"            = "#269d69",
                "Seurat + BPCells"  = "#1c8ca8",
                "rapids-singlecell" = "#7058be")
lib_order <- names(lib_colors)
method_to_group <- c(brisc = "brisc", scanpy = "Scanpy",
                     seurat = "Seurat + BPCells",
                     rapids = "rapids-singlecell")

op_order <- list(
  basic = c("Load data", "Quality control", "Normalization",
            "Feature selection", "PCA", "Nearest neighbors",
            "Clustering", "Embedding", "Find markers"),
  transfer = c("Load data", "Quality control", "Split data",
               "Normalization", "Feature selection", "PCA",
               "Transfer labels"),
  de = c("Load data", "Quality control", "Pseudobulk", "Filter",
         "Differential expression"),
  commands = c("Get expression by cell", "Get expression by gene",
               "Subset to one cell type",
               "Subset to highly variable genes",
               "Subsample to 10,000 cells", "Select categorical columns",
               "Split by cell type", "Concatenate cell types"))

op_wrap <- c(
  "Clustering"                      = "Cluster\n(5 resolutions)",
  "Get expression by cell"          = "Get expression\nby cell",
  "Get expression by gene"          = "Get expression\nby gene",
  "Subset to one cell type"         = "Subset to one\ncell type",
  "Subset to highly variable genes" = "Subset to highly\nvariable genes",
  "Subsample to 10,000 cells"       = "Subsample to\n10,000 cells",
  "Select categorical columns"      = "Select categorical\ncolumns",
  "Split by cell type"              = "Split by\ncell type",
  "Concatenate cell types"          = "Concatenate\ncell types")
op_lab <- setNames(unique(unlist(op_order)), unique(unlist(op_order)))
op_lab[names(op_wrap)] <- op_wrap

TIMEOUT_S <- 24 * 3600
TIMEOUT_ARROW <- 1.3
MIN_CLUSTER_SIZE <- 100

update_geom_defaults("text", list(family = font))
base_theme <- theme_bw(base_family = font) + theme(
  text = element_text(color = "black"),
  axis.text = element_text(color = "black"),
  plot.title = element_text(size = 11, hjust = 0.5),
  panel.grid = element_blank(),
  plot.margin = margin(6, 6, 6, 6))

fmt_time <- function(s) {
  case_when(s >= 3600 ~ paste0(round(s / 3600, 1), "h"),
            s >= 60   ~ paste0(round(s / 60, 1), "m"),
            s >= 0.1  ~ paste0(round(s, 1), "s"),
            .default  = paste0(round(s * 1000), "ms"))
}

make_label <- function(time_s, mult) {
  m <- if_else(mult < 10, as.character(round(mult, 1)),
               as.character(round(mult)))
  if_else(!is.na(mult) & is.finite(mult) & mult > 0,
          paste0(fmt_time(time_s), " (", m, "x)"), fmt_time(time_s))
}

make_title <- function(ds) {
  ggplot() + theme_void(base_family = font) +
    labs(title = ds_lab[ds]) +
    theme(plot.title = element_text(hjust = 0.5, size = 11,
                                    margin = margin(0, 0, 0, 0)),
          plot.margin = margin(0, 0, 0, 0))
}

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = width,
         height = height, dpi = 300, device = agg_png)
  ggsave(file.path(fig_dir, paste0(name, ".svg")), plot, width = width,
         height = height, device = svglite)
}

results <- list.files(file.path(work_dir, "output"), "_timer\\.csv$",
                      full.names = TRUE, recursive = TRUE) %>%
  keep(\(f) str_detect(f, "/rep\\d+/")) %>%
  map(\(f) fread(f) %>% mutate(.src = basename(f))) %>%
  rbindlist(fill = TRUE) %>% as_tibble() %>%
  mutate(gpu_run = str_detect(.src, "_gpu_timer\\.csv$"))

prepared <- results %>%
  filter(library != "brisc" |
           !str_detect(operation, "^Embedding \\(") |
           operation == "Embedding (LocalMAP)") %>%
  mutate(
    operation = if_else(operation == "Embedding (LocalMAP)", "Embedding",
                        operation),
    workflow = recode(test, manipulation = "commands") %>%
      factor(levels = c(names(wf_lab), "commands")),
    dataset = factor(dataset, datasets),
    group = case_when(
      library == "brisc" & gpu_run ~ "brisc (GPU node)",
      library == "brisc" & num_threads == -1 ~ "brisc (multi-threaded)",
      library == "brisc" ~ "brisc (single-threaded)",
      library == "scanpy" ~ "Scanpy",
      library == "rapids" ~ "rapids-singlecell",
      .default = "Seurat + BPCells") %>%
      factor(levels = c(grp_order, "rapids-singlecell",
                        "brisc (GPU node)", grp_gpu)),
    timed_out = FALSE)

impute_timeouts <- function(d, wf) {
  d_wf <- d %>% filter(workflow == wf)
  ops <- intersect(op_order[[wf]], unique(d_wf$operation))
  d_wf %>%
    distinct(group, dataset, rep) %>%
    crossing(operation = ops) %>%
    anti_join(d_wf %>% distinct(group, dataset, rep, operation),
              by = c("group", "dataset", "rep", "operation")) %>%
    slice_min(match(operation, ops), by = c(group, dataset, rep), n = 1) %>%
    left_join(d_wf %>%
                summarize(sum_dur = sum(duration), max_mem = max(memory),
                          .by = c(group, dataset, rep)),
              by = c("group", "dataset", "rep")) %>%
    transmute(group, dataset, rep, operation,
              duration = pmax(TIMEOUT_S - sum_dur, 0), memory = max_mem,
              workflow = factor(wf, levels = levels(d$workflow)),
              test = wf, timed_out = TRUE)
}

prepared <- prepared %>%
  bind_rows(map(names(wf_lab), \(wf) impute_timeouts(prepared, wf)) %>%
              list_rbind())

prepared <- prepared %>%
  anti_join(
    prepared %>%
      filter(workflow %in% names(wf_lab)) %>%
      summarize(incomplete = any(timed_out),
                .by = c(group, workflow, dataset, rep)) %>%
      filter(incomplete & any(!incomplete),
             .by = c(group, workflow, dataset)) %>%
      select(group, workflow, dataset, rep),
    by = c("group", "workflow", "dataset", "rep"))

prepared_main <- prepared %>%
  filter(group %in% grp_order) %>%
  mutate(group = factor(group, grp_order))

ref_total <- prepared_main %>%
  filter(workflow %in% names(wf_lab),
         group == "brisc (multi-threaded)") %>%
  summarize(total_s = sum(duration), .by = c(workflow, dataset, rep)) %>%
  summarize(ref = mean(total_s), .by = c(workflow, dataset))

totals <- prepared_main %>%
  filter(workflow %in% names(wf_lab)) %>%
  summarize(total_s = sum(duration), timed_out = any(timed_out),
            .by = c(group, workflow, dataset, rep)) %>%
  summarize(mean_s = mean(total_s), lo_s = min(total_s),
            hi_s = max(total_s), timed_out = any(timed_out), n = n(),
            .by = c(group, workflow, dataset)) %>%
  left_join(ref_total, by = c("workflow", "dataset")) %>%
  mutate(ref = if_else(group == "brisc (multi-threaded)", NA_real_, ref),
         div = ds_div[as.character(dataset)],
         dur_plot = mean_s / div,
         lo = if_else(n > 3, NA_real_, lo_s / div),
         hi = if_else(n > 3, NA_real_, hi_s / div),
         label = make_label(mean_s, mean_s / ref),
         label = if_else(timed_out, paste0(label, "\ntimed out"), label),
         label_x = if_else(is.na(hi), dur_plot, hi),
         label_x = if_else(timed_out,
                           pmax(label_x, dur_plot * TIMEOUT_ARROW),
                           label_x))

peak_mem <- prepared_main %>%
  filter(workflow %in% names(wf_lab)) %>%
  summarize(memory = max(memory), .by = c(group, workflow, dataset))

plot_bar <- function(wf, ds) {
  d <- totals %>% filter(workflow == wf, dataset == ds)
  pad <- max(d$dur_plot) * 0.04
  ggplot(d, aes(dur_plot, group, fill = group)) +
    geom_col(width = 0.7) +
    geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y",
                   linewidth = 0.3, colour = "grey25", na.rm = TRUE) +
    geom_segment(data = ~ filter(.x, timed_out),
                 aes(x = dur_plot * 0.98, xend = dur_plot * TIMEOUT_ARROW,
                     y = group, yend = group, color = group),
                 arrow = arrow(length = unit(0.06, "inches"),
                               type = "closed"),
                 linewidth = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(aes(x = label_x + pad, label = label), hjust = 0, size = 2) +
    scale_color_manual(values = grp_colors, guide = "none") +
    scale_fill_manual(values = grp_colors, guide = "none") +
    scale_x_continuous(
      breaks = scales::breaks_pretty(3),
      expand = expansion(mult = c(0, if (any(d$timed_out)) 0.9 else 0.5))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = if (wf == "de") paste0("Total runtime (", ds_unit[ds], ")"),
         y = if (ds == "SEAAD") str_replace(wf_lab[wf], " ", "\n")) +
    base_theme +
    theme(axis.title.x = element_text(size = 8, margin = margin(t = 5)),
          axis.title.y = element_text(size = 11, margin = margin(r = 10)),
          axis.text.x = element_text(size = 8),
          axis.text.y = if (ds == "SEAAD") element_text(size = 8)
                        else element_blank(),
          axis.ticks.y = if (ds == "SEAAD") element_line()
                         else element_blank(),
          plot.margin = margin(0, 0, 2, if (ds == "SEAAD") 4 else 8))
}

plot_mem <- function(wf, ds) {
  d <- peak_mem %>% filter(workflow == wf, dataset == ds)
  ggplot(d, aes(memory, group)) +
    geom_col(fill = "#d9d9d9", width = 0.7) +
    geom_text(aes(x = 0, label = paste0(round(memory), " GB")),
              hjust = 0, nudge_x = max(d$memory) * 0.08, size = 2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.5))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = if (wf == "de") "Peak memory", y = NULL) +
    base_theme +
    theme(axis.title.x = element_text(size = 8, margin = margin(t = 5)),
          axis.text = element_blank(), axis.ticks = element_blank(),
          plot.margin = margin(0, if (ds == "PanSci") 4 else 0, 2, -1))
}

total_panels <- c(
  map(datasets, make_title),
  map(names(wf_lab), \(wf) {
    map(datasets, \(ds) list(plot_bar(wf, ds), plot_mem(wf, ds))) %>%
      list_flatten()
  }) %>% list_flatten())
names(total_panels) <- LETTERS[seq_along(total_panels)]

save_fig("fig_total_runtime",
         wrap_plots(total_panels, design = "AABBCC\nDEFGHI\nJKLMNO\nPQRSTU",
                    heights = c(0.15, 10, 10, 10),
                    widths = rep(c(4, 1.3), 3)),
         width = 11, height = 4.4)

compute_steps <- function(d, ref_group, annotate = "others") {
  by_op <- c("operation", "workflow", "dataset")
  agg <- d %>%
    summarize(dur_min = min(duration), dur_max = max(duration),
              duration = mean(duration), n = n(), memory = max(memory),
              timed_out = any(timed_out),
              .by = c(group, operation, workflow, dataset)) %>%
    left_join(d %>%
                filter(group == ref_group) %>%
                summarize(ref = mean(duration), .by = all_of(by_op)),
              by = by_op) %>%
    left_join(d %>%
                filter(group != ref_group) %>%
                summarize(other = mean(duration), .by = all_of(by_op)),
              by = by_op)
  agg %>%
    mutate(mult = if (annotate == "ref")
                    if_else(group == ref_group, other / duration, NA_real_)
                  else if_else(group == ref_group, NA_real_, duration / ref),
           div = ds_div[as.character(dataset)],
           dur_plot = duration / div,
           label = make_label(duration, mult),
           label = if_else(timed_out, paste0(label, "\ntimed out"),
                           label)) %>%
    mutate(dur_plot = pmax(dur_plot, max(dur_plot) * 0.02),
           .by = c(operation, workflow, dataset)) %>%
    mutate(lo = if_else(n > 3, NA_real_,
                        pmax(dur_plot - (duration - dur_min) / div, 0)),
           hi = if_else(n > 3, NA_real_,
                        dur_plot + (dur_max - duration) / div),
           label_x = if_else(is.na(hi), dur_plot, hi),
           label_x = if_else(timed_out,
                             pmax(label_x, dur_plot * TIMEOUT_ARROW),
                             label_x)) %>%
    mutate(label_x = label_x + max(dur_plot) * 0.05,
           .by = c(operation, workflow, dataset))
}

steps_main <- compute_steps(prepared_main, "brisc (multi-threaded)")

step_panel <- function(d, ds, group_order, legend_title, is_first,
                       right_margin) {
  ggplot(d %>% filter(dataset == ds), aes(dur_plot, group, fill = group)) +
    geom_col(width = 0.8) +
    geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y",
                   linewidth = 0.25, colour = "grey25", na.rm = TRUE) +
    geom_segment(data = ~ filter(.x, timed_out),
                 aes(x = dur_plot * 0.98, xend = dur_plot * TIMEOUT_ARROW,
                     y = group, yend = group, color = group),
                 arrow = arrow(length = unit(0.06, "inches"),
                               type = "closed"),
                 linewidth = 0.8, inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(aes(x = label_x, label = label), hjust = 0, size = 2) +
    facet_wrap(~ operation, ncol = 1, scales = "free_x",
               strip.position = "left") +
    scale_color_manual(values = grp_colors, guide = "none") +
    scale_fill_manual(values = grp_colors, limits = group_order,
                      name = legend_title,
                      guide = if (is_first) guide_legend(reverse = TRUE)
                              else "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 1.1))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Per-step runtime", y = NULL) +
    base_theme +
    theme(strip.text.y.left = if (is_first)
            element_text(angle = 0, hjust = 1, size = 8)
            else element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_text(size = 8, margin = margin(t = 5)),
          panel.spacing.y = unit(-0.5, "pt"),
          plot.margin = margin(0, right_margin, 2,
                               if (is_first) 4 else 8))
}

mem_panel <- function(d, ds, is_last) {
  d_ds <- d %>% filter(dataset == ds)
  mem_max <- max(d_ds$memory)
  ggplot(d_ds, aes(memory, group)) +
    geom_col(fill = "#d9d9d9", width = 0.8) +
    geom_text(aes(x = 0, label = paste0(round(memory), " GB")), hjust = 0,
              nudge_x = mem_max * 0.08, size = 2) +
    facet_wrap(~ operation, ncol = 1, scales = "free_x",
               strip.position = "right") +
    scale_x_continuous(limits = c(0, mem_max * 1.5),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Peak\nmemory", y = NULL) +
    base_theme +
    theme(strip.text = element_blank(), strip.background = element_blank(),
          strip.placement = "outside",
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.ticks.length.y = unit(0, "pt"),
          axis.title.x = element_text(size = 8, margin = margin(t = 5)),
          panel.spacing.y = unit(-0.5, "pt"),
          plot.margin = margin(0, if (is_last) 4 else 0, 2, -1))
}

build_step_fig <- function(wf, steps_df, group_order, legend_title,
                           memory = TRUE, titles = TRUE) {
  ops <- intersect(op_order[[wf]],
                   unique(steps_df$operation[steps_df$workflow == wf]))
  d <- steps_df %>%
    filter(workflow == wf, operation %in% ops) %>%
    mutate(operation = factor(operation, ops, labels = op_lab[ops]))

  title_row <- if (titles) map(datasets, make_title) else list()
  steps <- imap(datasets, \(ds, j)
    step_panel(d, ds, group_order, legend_title, j == 1,
               if (memory) 0 else if (j == length(datasets)) 4 else 8))

  if (memory) {
    panels <- c(title_row, list_flatten(map(seq_along(datasets), \(j)
      list(steps[[j]], mem_panel(d, datasets[j], j == length(datasets))))))
    design <- if (titles) "AABBCC\nDEFGHI" else "ABCDEF"
    widths <- rep(c(4, 1.3), 3)
  } else {
    panels <- c(title_row, steps)
    design <- if (titles) "ABC\nDEF" else "ABC"
    widths <- rep(4, 3)
  }
  names(panels) <- LETTERS[seq_along(panels)]

  wrap_plots(panels, design = design,
             heights = if (titles) c(0.15, 10) else 10,
             widths = widths) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9))
}

walk(names(op_order), \(wf) {
  n <- length(intersect(op_order[[wf]],
                        unique(steps_main$operation[
                          steps_main$workflow == wf])))
  save_fig(paste0("fig_", wf, "_steps"),
           build_step_fig(wf, steps_main, grp_order,
                          "192 CPUs, 749 GB RAM"),
           width = 11, height = max(3, n * 0.85))
})

knn_data <- fread(file.path(work_dir, "output",
                            "knn_recall_summary.csv")) %>%
  as_tibble() %>%
  mutate(dataset = factor(dataset, datasets),
         lib = factor(method_to_group[method], levels = lib_order)) %>%
  inner_join(
    prepared %>%
      filter(workflow == "basic", operation == "Nearest neighbors",
             !group %in% c("brisc (single-threaded)",
                           "brisc (GPU node)")) %>%
      summarize(nn_s = mean(duration), .by = c(group, dataset)) %>%
      transmute(lib = factor(str_remove(as.character(group),
                                        " \\(multi-threaded\\)"),
                             lib_order),
                dataset, nn_s,
                hw = if_else(group == "rapids-singlecell", "GPU node",
                             "192-CPU node")),
    by = c("lib", "dataset"))

knn_k <- paste(sort(unique(knn_data$k)), collapse = ", ")
key_hw <- coalesce(deframe(distinct(knn_data, lib, hw))[lib_order],
                   "192-CPU node")

plot_knn <- knn_data %>%
  select(lib, dataset, nn_s, hw, mean, cvar25) %>%
  pivot_longer(c(mean, cvar25), names_to = "stat", values_to = "recall") %>%
  mutate(stat = factor(stat, c("mean", "cvar25"),
                       labels = c("All cells", "Bottom 25% of cells"))) %>%
  ggplot(aes(nn_s, recall)) +
  geom_segment(data = knn_data, inherit.aes = FALSE,
               aes(x = nn_s, xend = nn_s, y = cvar25, yend = mean,
                   colour = lib),
               linewidth = 0.7, alpha = 0.55) +
  geom_point(aes(colour = lib, alpha = stat, shape = hw), size = 2.8) +
  ggrepel::geom_text_repel(
    data = knn_data, inherit.aes = FALSE,
    aes(nn_s * 1.3, mean, label = sprintf("%.4f", mean), colour = lib),
    size = 2.3, hjust = 0, direction = "y", box.padding = 0.12,
    min.segment.length = Inf, seed = 1, show.legend = FALSE) +
  facet_wrap(~ dataset, nrow = 1, scales = "free",
             labeller = as_labeller(ds_lab)) +
  scale_colour_manual(
    values = lib_colors, limits = lib_order, name = NULL,
    labels = if_else(key_hw == "GPU node", paste0(lib_order, " (GPU)"),
                     recode(lib_order, brisc = "brisc (multi-threaded)")),
    guide = guide_legend(order = 1, override.aes = list(
      alpha = 1, shape = if_else(key_hw == "GPU node", 17, 16),
      linetype = 0))) +
  scale_shape_manual(values = c("192-CPU node" = 16, "GPU node" = 17),
                     guide = "none") +
  scale_alpha_manual(
    values = c("All cells" = 1, "Bottom 25% of cells" = 0.4), name = NULL,
    guide = guide_legend(order = 2, override.aes = list(
      colour = "grey25", size = 2.8, shape = 16))) +
  scale_x_log10(breaks = c(1, 10, 60, 600, 3600, 36000),
                labels = fmt_time,
                expand = expansion(mult = c(0.06, 0.2))) +
  scale_y_continuous(breaks = scales::breaks_pretty(4),
                     expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = "Nearest-neighbor step runtime",
       y = paste0("Mean kNN recall at k = ", knn_k, "\n(vs exact kNN)")) +
  base_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 9)),
        axis.title.y = element_text(margin = margin(r = 9)),
        legend.position = "right",
        legend.spacing.y = unit(2, "pt"),
        panel.grid.major = element_line(color = "grey93"))

save_fig("fig_knn_recall", plot_knn, width = 11.5, height = 3.4)

gpu_rename <- c("rapids-singlecell" = "rapids-singlecell (GPU)",
                "brisc (GPU node)"  = "brisc (CPU)")
steps_gpu <- prepared %>%
  filter(group %in% names(gpu_rename)) %>%
  mutate(group = factor(gpu_rename[as.character(group)], levels = grp_gpu)) %>%
  compute_steps("brisc (CPU)", annotate = "ref")

gpu_totals <- prepared %>%
  filter(group %in% names(gpu_rename), workflow == "basic") %>%
  mutate(group = factor(gpu_rename[as.character(group)], levels = grp_gpu)) %>%
  summarize(total_s = sum(duration), .by = c(group, dataset, rep)) %>%
  summarize(mean_s = mean(total_s), lo_s = min(total_s), hi_s = max(total_s),
            .by = c(group, dataset)) %>%
  mutate(mult = if_else(group == "brisc (CPU)",
                        mean_s[group != "brisc (CPU)"] / mean_s, NA_real_),
         .by = dataset) %>%
  mutate(label = make_label(mean_s, mult))

gpu_total_panel <- function(ds) {
  d <- gpu_totals %>% filter(dataset == ds)
  unit <- if (max(d$mean_s) >= 3600) "h" else
          if (max(d$mean_s) >= 600) "m" else "s"
  div <- c(h = 3600, m = 60, s = 1)[[unit]]
  d <- d %>% mutate(dur_plot = mean_s / div, lo = lo_s / div,
                    hi = hi_s / div)
  ggplot(d, aes(dur_plot, group, fill = group)) +
    geom_col(width = 0.7) +
    geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y",
                   linewidth = 0.3, colour = "grey25") +
    geom_text(aes(x = hi + max(d$dur_plot) * 0.04, label = label),
              hjust = 0, size = 2) +
    scale_fill_manual(values = grp_colors, guide = "none") +
    scale_x_continuous(breaks = scales::breaks_pretty(3),
                       expand = expansion(mult = c(0, 0.5))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = paste0("Total runtime (", unit, ")"), y = NULL) +
    base_theme +
    theme(axis.title.x = element_text(size = 8, margin = margin(t = 5)),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(0, if (ds == "PanSci") 4 else 8, 2,
                               if (ds == "SEAAD") 4 else 8))
}

gpu_panels <- c(map(datasets, make_title), map(datasets, gpu_total_panel))
names(gpu_panels) <- LETTERS[seq_along(gpu_panels)]

n_gpu <- length(intersect(op_order$basic,
                          unique(steps_gpu$operation[
                            steps_gpu$workflow == "basic"])))
gpu_panels$A <- gpu_panels$A + labs(tag = "a")
gpu_steps_fig <- build_step_fig("basic", steps_gpu, grp_gpu,
                                "96 CPUs,\n4× H100 GPUs,\n752 GB RAM",
                                memory = FALSE, titles = FALSE)
gpu_steps_fig[[1]] <- gpu_steps_fig[[1]] + labs(tag = "b")

save_fig("fig_basic_gpu_steps",
         wrap_plots(gpu_panels, design = "AABBCC\nDDEEFF",
                    heights = c(0.15, 4)) / gpu_steps_fig +
           plot_layout(heights = c(1.5, n_gpu * 0.42)) &
           theme(plot.tag = element_text(size = 12, face = "bold"),
                 plot.tag.position = c(0, 1)),
         width = 10.5, height = max(2, n_gpu * 0.42) + 1.7)

ari_data <- fread(file.path(work_dir, "output",
                            "cluster_ari_summary.csv")) %>%
  as_tibble() %>%
  rename(k = all_of(str_subset(names(.), "^n_clusters_ge_\\d+$"))) %>%
  filter(label == "cell_type") %>%
  mutate(dataset = factor(dataset, datasets),
         group = factor(method_to_group[method], levels = lib_order))

plot_ari <- ggplot(ari_data, aes(k, balanced_ari, color = group)) +
  geom_vline(data = ari_data %>%
               summarize(n = max(n_classes), .by = dataset),
             inherit.aes = FALSE, aes(xintercept = n),
             linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.6) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x",
             labeller = as_labeller(ds_lab)) +
  scale_color_manual(values = lib_colors, limits = lib_order, name = NULL) +
  scale_x_log10(breaks = c(10, 15, 25, 40, 60, 100, 150)) +
  scale_y_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Number of clusters",
       y = "Balanced ARI\n(compared to cell type labels)") +
  base_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 9)),
        axis.title.y = element_text(margin = margin(r = 9)),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey93"))

save_fig("fig_cluster_ari", plot_ari, width = 11, height = 3.4)

embed_paths <- c(brisc  = "basic_brisc_%s_-1_embedding.csv",
                 scanpy = "basic_scanpy_%s_embedding.csv",
                 seurat = "basic_seurat_%s_embedding.csv",
                 rapids = "basic_rapids_%s_gpu_embedding.csv")

panel_spec <- tibble(
  lib = c(rep("brisc", 4), "scanpy", "seurat", "rapids"),
  method = c("LocalMAP", "PaCMAP", "UMAP", "UMAP hogwild", rep("UMAP", 3)),
  col1 = c("localmap_1", "pacmap_1", "umap_1", "umap_hogwild_1",
           rep("embed_1", 3)),
  col2 = str_replace(col1, "_1$", "_2"),
  op = c("Embedding (LocalMAP)", "Embedding (PaCMAP)", "Embedding (UMAP)",
         "Embedding (UMAP hogwild)", rep("Embedding", 3)))

embed_time <- function(lib_val, op_val, ds_val) {
  results %>%
    filter(test == "basic", dataset == ds_val, rep == 1,
           library == lib_val, operation == op_val,
           !(library == "brisc" & num_threads != -1),
           !(library == "brisc" & str_detect(.src, "_gpu_timer"))) %>%
    pull(duration) %>% first()
}

build_embed_row <- function(ds) {
  files <- file.path(work_dir, "output", "rep1", sprintf(embed_paths, ds))
  names(files) <- names(embed_paths)

  set.seed(1)
  raw <- imap(files, \(f, lib_name) {
    cols <- panel_spec %>% filter(lib == lib_name) %>% select(col1, col2)
    df <- fread(f, select = c("cluster_res_1.0", "cell_type_broad",
                              unique(unlist(cols)))) %>% as_tibble()
    df[sample.int(nrow(df)), ]
  })

  types <- sort(unique(unlist(map(raw, \(d) unique(d$cell_type_broad)))))
  hues <- setNames(head(seq(0, 360, length.out = length(types) + 1), -1),
                   types)

  pal <- imap(raw, \(d, lib_name) {
    d %>%
      count(cluster_res_1.0, cell_type_broad) %>%
      mutate(n_cells = sum(n), .by = cluster_res_1.0) %>%
      slice_max(n, by = cluster_res_1.0, with_ties = FALSE) %>%
      arrange(as.numeric(cluster_res_1.0)) %>%
      mutate(big = n_cells >= MIN_CLUSTER_SIZE) %>%
      mutate(rank = cumsum(big), n_cl = sum(big), .by = cell_type_broad) %>%
      mutate(l = if_else(big, 85 - 40 * (rank - 1) / pmax(n_cl - 1, 1), 65),
             color = hcl(hues[cell_type_broad], 55, pmin(pmax(l, 45), 85)),
             key = paste0(lib_name, "_", cluster_res_1.0))
  })
  keys <- unlist(map(pal, "key"), use.names = FALSE)
  colors <- setNames(unlist(map(pal, "color"), use.names = FALSE),
                     seq_along(keys))

  panel_info <- panel_spec %>%
    mutate(panel_key = paste0(lib, "_", method),
           duration = map2_dbl(lib, op, embed_time, ds),
           label = paste0(method_to_group[lib], " (", method, ", ",
                          fmt_time(duration), ")"))

  embed_data <- panel_spec %>%
    pmap(\(lib, method, col1, col2, op) {
      d <- raw[[lib]]
      tibble(embed_1 = d[[col1]], embed_2 = d[[col2]],
             cluster_id = factor(match(paste0(lib, "_", d$cluster_res_1.0),
                                       keys), seq_along(keys)),
             panel_label = factor(paste0(lib, "_", method),
                                  panel_info$panel_key, panel_info$label))
    }) %>%
    bind_rows()

  legend_data <- imap(raw, \(d, lib_name) distinct(d, cell_type_broad)) %>%
    bind_rows() %>%
    distinct(cell_type_broad) %>%
    mutate(embed_1 = NA_real_, embed_2 = NA_real_,
           panel_label = factor(panel_info$label[1], panel_info$label))

  ggplot(embed_data, aes(embed_1, embed_2)) +
    geom_scattermore(aes(color = cluster_id), pointsize = 1,
                     pixels = c(700, 700)) +
    geom_point(data = legend_data, aes(fill = cell_type_broad), shape = 22,
               size = 0, stroke = 0, na.rm = TRUE) +
    scale_color_manual(values = colors, guide = "none") +
    scale_fill_manual(values = setNames(hcl(hues[types], 55, 60), types),
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = 3))) +
    facet_wrap(~ panel_label, nrow = 1, scales = "free") +
    labs(title = ds_lab[ds], x = NULL, y = NULL) +
    base_theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 8),
          legend.key.spacing.y = unit(1, "pt"))
}

save_fig("fig_embeddings",
         wrap_plots(map(datasets, build_embed_row), ncol = 1),
         width = 15, height = 2.6 * length(datasets))

acc_data <- list.files(file.path(work_dir, "output"),
                       "^transfer_.*_accuracy\\.csv$", full.names = TRUE,
                       recursive = TRUE) %>%
  keep(\(f) str_detect(f, "/rep\\d+/")) %>%
  map(\(f) fread(f) %>% as_tibble() %>% mutate(.src = basename(f))) %>%
  bind_rows() %>%
  filter(cell_type != "Total") %>%
  mutate(library = str_extract(.src, "(?<=^transfer_)[^_]+"),
         nt = suppressWarnings(as.integer(
           str_extract(.src, "(?<=_)(-?\\d+)(?=_accuracy)")))) %>%
  filter(library != "brisc" | nt == -1) %>%
  mutate(dataset = factor(str_extract(.src, "SEAAD|Parse|PanSci"), datasets),
         lib = factor(method_to_group[library], levels = lib_order))

acc_type <- acc_data %>%
  summarize(pct = mean(percent_correct), .by = c(lib, dataset, cell_type))

acc_note <- acc_data %>%
  summarize(pct = sum(n_correct) / sum(n_total) * 100,
            .by = c(lib, dataset)) %>%
  transmute(lib, dataset, text = sprintf("%s: %.1f%%", lib, pct)) %>%
  bind_rows(
    prepared %>%
      filter(workflow == "transfer") %>%
      mutate(lib = factor(if_else(str_detect(group, "^brisc"), "brisc",
                                  as.character(group)), lib_order)) %>%
      summarize(timed_out = any(timed_out), .by = c(lib, dataset)) %>%
      filter(timed_out) %>%
      transmute(lib, dataset, text = paste0(lib, ": timed out"))) %>%
  arrange(dataset, lib) %>%
  mutate(y = 0.23 - 0.07 * (row_number() - 1), .by = dataset)

acc_head <- acc_note %>%
  summarize(y_head = 0.30, y_count = min(y) - 0.07, .by = dataset) %>%
  left_join(acc_type %>%
              summarize(n_types = n_distinct(cell_type), .by = dataset),
            by = "dataset")

plot_acc <- ggplot(acc_type, aes(pct, colour = lib)) +
  stat_ecdf(aes(y = after_stat(1 - ecdf)), geom = "step",
            linewidth = 0.6) +
  geom_text(data = acc_note, inherit.aes = FALSE,
            aes(x = 2, y = y, label = text, colour = lib), hjust = 0,
            size = 2.3, show.legend = FALSE) +
  geom_text(data = acc_head, inherit.aes = FALSE, aes(x = 2, y = y_head),
            label = "Accuracy over all cells", hjust = 0, size = 2.3,
            colour = "grey30") +
  geom_text(data = acc_head, inherit.aes = FALSE,
            aes(x = 2, y = y_count, label = paste(n_types, "cell types")),
            hjust = 0, size = 2.3, colour = "grey30") +
  facet_wrap(~ dataset, nrow = 1, labeller = as_labeller(ds_lab)) +
  scale_colour_manual(values = lib_colors, limits = lib_order,
                      breaks = setdiff(lib_order, "rapids-singlecell"),
                      name = NULL) +
  scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Accuracy per cell type (%)",
       y = "Fraction of cell types\nat least this accurate") +
  base_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 9)),
        axis.title.y = element_text(margin = margin(r = 9)),
        legend.position = "right",
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "grey93"))

save_fig("fig_transfer_accuracy", plot_acc, width = 11, height = 3.4)

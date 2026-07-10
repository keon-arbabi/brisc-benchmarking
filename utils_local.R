suppressPackageStartupMessages({
  library(dplyr)
  library(processx)
})

.this_file_path <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
.MONITOR_MEM_SH_PATH <- file.path(.this_file_path, "monitor_mem.sh")

POLLING_INTERVAL <- '0.05'

.TIME_CONVERSIONS <- c(
  "s" = 1, "ms" = 1000, "us" = 1e6, "\u00b5s" = 1e6,
  "ns" = 1e9, "m" = 1/60, "h" = 1/3600, "d" = 1/86400
)

MemoryTimer = function(silent = TRUE, csv_path = NULL, csv_columns = NULL,
                       unit = "s", summary_unit = NULL) {
  env = environment()
  env$timings = list()
  pid = Sys.getpid()
  env$.shutdown_done = FALSE
  env$.csv_path = csv_path
  env$.csv_columns = csv_columns
  env$.unit = unit
  env$.summary_unit = if (!is.null(summary_unit)) summary_unit else unit

  with_timer = function(message, expr, exclude = FALSE) {
    if (!silent) cat(paste0(message, '...\n'))
    result = NULL
    aborted = FALSE
    # Write monitor output to a temp file instead of a pipe. R is single-
    # threaded: during long-running expressions (e.g. rep() for 150+ GiB),
    # the 64 KB pipe buffer fills and the monitor blocks on write, missing
    # the peak RSS. A file has no such limit.
    tmpf <- tempfile(pattern = "memmon_", fileext = ".csv")
    curr_process <- process$new(
      command = .MONITOR_MEM_SH_PATH,
      args = c("-p", as.character(pid), "-i", POLLING_INTERVAL),
      stdout = tmpf
    )
    # Sync: wait for first sample before timing (file must have content)
    deadline <- proc.time()[['elapsed']] + 5
    repeat {
      if (file.exists(tmpf) && file.info(tmpf)$size > 0) break
      if (proc.time()[['elapsed']] > deadline) break
      Sys.sleep(0.005)
    }
    start = Sys.time()

    # Store error for re-throwing
    error_obj <- NULL

    tryCatch({
      result = invisible(eval(substitute(expr), parent.frame()))
    }, interrupt = function(i) {
      aborted <<- TRUE
      error_obj <<- i
    }, error = function(e) {
      aborted <<- TRUE
      error_obj <<- e
    }, finally = {
      duration = as.numeric(difftime(Sys.time(), start, units = 'secs'))
      curr_process$signal(15L)  # SIGTERM
      curr_process$wait(timeout = 5000)
      if (curr_process$is_alive()) curr_process$kill()
      stdout_output <- tryCatch(
        paste(readLines(tmpf, warn = FALSE), collapse = "\n"),
        error = function(e) "")
      unlink(tmpf)
      con <- textConnection(stdout_output)
      # Handle empty output
      df <- tryCatch({
        read.csv(con, header = FALSE, sep = ",", strip.white = TRUE,
          col.names = c("Integer", "Percentage"))
      }, error = function(e) {
        data.frame(Integer = 0, Percentage = 0.0)
      })
      close(con)
      # Ensure valid data
      if (nrow(df) == 0 || all(is.na(df$Integer))) {
        peak_mem <- 0
        percent <- 0.0
      } else {
        peak_mem <- max(df$Integer, na.rm = TRUE)
        percent <- max(df$Percentage, na.rm = TRUE)
        if (is.infinite(peak_mem)) peak_mem <- 0
        if (is.infinite(percent)) percent <- 0.0
      }

      new_memory <- round(peak_mem / 1024 / 1024, 2)
      new_percent <- round(percent, 2)

      if (message %in% names(env$timings)) {
        env$timings[[message]]$duration <-
          env$timings[[message]]$duration + duration
        env$timings[[message]]$max_mem <-
          max(env$timings[[message]]$max_mem, new_memory)
        env$timings[[message]]$mem_percent <-
          max(env$timings[[message]]$mem_percent, new_percent)
        env$timings[[message]]$aborted <-
          env$timings[[message]]$aborted || aborted
      } else {
        env$timings[[message]] = list(
          duration = duration,
          max_mem = new_memory,
          mem_percent = new_percent,
          aborted = aborted,
          exclude = exclude
        )
      }

      if (!silent) {
        time_str = format_time(duration)
        status = if (aborted) 'aborted after' else 'took'
        cat(sprintf(
          '%s %s %s using %.2f GiB\n\n',
          message, status, time_str, new_memory))
      }
      write_csv_now()
      gc()
    })
    # Re-throw after cleanup, preserving the original condition class.
    # stop() only accepts error conditions; for interrupts use the
    # base calling handler approach.
    if (!is.null(error_obj)) {
      if (inherits(error_obj, "error")) {
        stop(error_obj)
      } else {
        # Interrupt (e.g. SLURM TERM->INT): the finally block already flushed
        # the timer CSV. signalCondition() is a silent no-op with no outer
        # handler, which would let later steps run on partial state, so halt
        # the whole job with a non-zero status (matching Python's behavior).
        quit(save = "no", status = 1L, runLast = TRUE)
      }
    }

    return(invisible(result))
  }

  format_time = function(duration, unit = NULL) {
    duration = as.numeric(duration)

    if (!is.null(unit)) {
      if (!(unit %in% names(.TIME_CONVERSIONS)))
        stop("Unsupported unit: ", unit)
      converted = duration * .TIME_CONVERSIONS[[unit]]
      return(paste0(format(converted, scientific = FALSE), unit))
    }

    units = list(
      list(threshold = 86400, suffix = 'd'),
      list(threshold = 3600, suffix = 'h'),
      list(threshold = 60, suffix = 'm'),
      list(threshold = 1, suffix = 's'),
      list(threshold = 0.001, suffix = 'ms'),
      list(threshold = 0.000001, suffix = 'µs'),
      list(threshold = 0.000000001, suffix = 'ns')
    )

    parts = c()
    for (unit in units) {
      threshold = unit$threshold
      suffix = unit$suffix
      if (duration >= threshold ||
          (length(parts) == 0 && threshold == 0.000000001)) {
        if (threshold >= 1) {
          value = as.integer(duration %/% threshold)
          duration = duration %% threshold
        } else {
          value = as.integer((duration / threshold) %% 1000)
        }
        if (value > 0 ||
            (length(parts) == 0 && threshold == 0.000000001)) {
          parts = c(parts, paste0(value, suffix))
        }
        if (length(parts) == 2) break
      }
    }

    if (length(parts) > 0) {
      paste(parts, collapse = ' ')
    } else {
      'less than 1ns'
    }
  }

  print_summary = function(sort = FALSE, unit = NULL) {
    cat('\n--- Timing Summary ---\n')

    if (length(env$timings) == 0) {
      cat('no timings recorded.\n')
      return(invisible(NULL))
    }

    if (sort) {
      durations = sapply(env$timings, function(x) x$duration)
      items = names(env$timings)[order(durations, decreasing = TRUE)]
    } else {
      items = names(env$timings)
    }

    total_time = sum(sapply(env$timings, function(x)
      if (isTRUE(x$exclude)) 0 else x$duration))

    # Create table headers
    duration_header = if (!is.null(unit)) {
      paste0('Duration (', unit, ')')
    } else {
      'Duration'
    }
    headers = c('Operation', 'Status', duration_header, '% of Total',
                'Memory (GiB)', '% of Avail')

    # Create table data
    table_data = list()
    for (i in seq_along(items)) {
      msg = items[i]
      info = env$timings[[msg]]
      duration = info$duration
      percentage = if (total_time > 0) {
        (duration / total_time) * 100
      } else {
        0
      }
      max_mem = info$max_mem
      mem_percent = info$mem_percent
      status = if (info$aborted) 'aborted' else 'completed'
      time_str = format_time(duration, unit)

      excluded <- isTRUE(info$exclude)
      table_data[[i]] = c(
        msg,
        status,
        time_str,
        if (excluded) 'N/A' else sprintf('%.2f%%', percentage),
        sprintf('%.2f', max_mem),
        if (excluded) 'N/A' else sprintf('%.2f%%', mem_percent)
      )
    }

    cat(sprintf('%-30s %-10s %-15s %-10s %-13s %-10s\n',
                headers[1], headers[2], headers[3], headers[4],
                headers[5], headers[6]))
    cat(paste(rep('-', 90), collapse = ''), '\n')

    for (row in table_data) {
      cat(sprintf('%-30s %-10s %-15s %-10s %-13s %-10s\n',
                  row[1], row[2], row[3], row[4], row[5], row[6]))
    }

    cat(sprintf('\nTotal time: %s\n', format_time(total_time, unit)))
  }

  to_dataframe = function(sort = TRUE, unit = NULL) {
    if (length(env$timings) == 0) {
      return(data.frame(
        operation = character(0), duration = numeric(0),
        duration_unit = character(0), aborted = logical(0),
        percentage = numeric(0), memory = numeric(0),
        memory_unit = character(0), percent_mem = numeric(0)
      ))
    }

    if (sort) {
      durations = sapply(env$timings, function(x) x$duration)
      items = names(env$timings)[order(durations, decreasing = TRUE)]
    } else {
      items = names(env$timings)
    }

    total = sum(sapply(env$timings, function(x)
      if (isTRUE(x$exclude)) 0 else x$duration))

    ops = items
    durs = sapply(items, function(msg) env$timings[[msg]]$duration)

    if (!is.null(unit)) {
      if (!(unit %in% names(.TIME_CONVERSIONS)))
        stop("Unsupported unit: ", unit)
      durs = durs * .TIME_CONVERSIONS[[unit]]
      duration_unit = unit
    } else {
      duration_unit = "s"
    }

    aborts = sapply(items, function(msg) env$timings[[msg]]$aborted)
    pcts = sapply(items, function(msg) {
      if (isTRUE(env$timings[[msg]]$exclude)) NA_real_
      else if (total > 0) (env$timings[[msg]]$duration / total) * 100
      else 0
    })
    memory = sapply(items, function(msg) env$timings[[msg]]$max_mem)
    memory_unit = rep("GiB", length(items))
    percent_mem = sapply(items, function(msg) {
      if (isTRUE(env$timings[[msg]]$exclude)) NA_real_
      else env$timings[[msg]]$mem_percent
    })

    data.frame(
      operation = ops,
      duration = durs,
      duration_unit = duration_unit,
      aborted = aborts,
      percentage = pcts,
      memory = memory,
      memory_unit = memory_unit,
      percent_mem = percent_mem
    )
  }

  write_csv_now = function() {
    if (is.null(env$.csv_path) || length(env$timings) == 0)
      return(invisible(NULL))
    df <- to_dataframe(sort = FALSE, unit = env$.unit)
    if (!is.null(env$.csv_columns)) {
      for (col_name in names(env$.csv_columns)) {
        df[[col_name]] <- env$.csv_columns[[col_name]]
      }
    }
    write.csv(df, env$.csv_path, row.names = FALSE)
  }

  shutdown = function() {
    if (env$.shutdown_done) return(invisible(NULL))
    env$.shutdown_done <- TRUE
    if (length(env$timings) == 0) return(invisible(NULL))
    print_summary(unit = env$.summary_unit)
    write_csv_now()
  }

  # Auto-flush timings on R exit (e.g. SLURM timeout)
  shutdown_ref <- new.env(parent = emptyenv())
  shutdown_ref$fn <- shutdown
  reg.finalizer(shutdown_ref, function(e) e$fn(), onexit = TRUE)

  structure(list(
    with_timer = with_timer,
    print_summary = print_summary,
    to_dataframe = to_dataframe,
    shutdown = shutdown
  ), class = "TimerCollection")
}

system_info <- function() {
  hostname <- Sys.info()['nodename']
  user <- Sys.getenv("USER")
  if (user == "") user <- "N/A"

  # Get CPU cores
  cpu_task <- Sys.getenv("SLURM_CPUS_PER_TASK")
  cpu_cores <- Sys.getenv("SLURM_CPUS_ON_NODE")

  if (nchar(cpu_task) > 0) {
    cpu_cores <- cpu_task
  } else if (nchar(cpu_cores) == 0) {
    tryCatch({
      if (requireNamespace("parallel", quietly = TRUE)) {
        cpu_cores <- parallel::detectCores()
      } else {
        cpu_cores <- "N/A"
      }
    }, error = function(e) {
      cpu_cores <<- "N/A"
    })
  }

  # Get Memory
  mem_gb_str <- "N/A"
  mem_mb_str <- Sys.getenv("SLURM_MEM_PER_NODE")

  if (nchar(mem_mb_str) > 0) {
    # First, try SLURM environment variable
    mem_mb <- suppressWarnings(as.numeric(mem_mb_str))
    if (!is.na(mem_mb)) {
      mem_gb_str <- sprintf("%.1f GB", mem_mb / 1024)
    }
  } else {
    # As a fallback, read from /proc/meminfo
    tryCatch({
      meminfo_lines <- readLines("/proc/meminfo")
      for (line in meminfo_lines) {
        if (grepl("^MemTotal:", line)) {
          mem_kb <- as.numeric(strsplit(line, "\\s+")[[1]][2])
          mem_gb_str <- sprintf("%.1f GB", mem_kb / 1024 / 1024)
          break
        }
      }
    }, error = function(e) {
      # Keep default "N/A"
    })
  }

  cat("\n--- User Resource Allocation ---\n")
  cat(sprintf("Node: %s\n", hostname))
  cat(sprintf("User: %s\n", user))
  cat(sprintf("CPU Cores Allocated: %s\n", cpu_cores))
  cat(sprintf("Memory Allocated: %s\n", mem_gb_str))
  cat(sprintf("R.version=%s\n", R.version.string))
  cat("\n")
}

read_h5ad_obs <- function(path) {
  # Load hdf5r only when needed
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' required for reading h5ad files")
  }

  h5_file <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5_file$close_all(), add = TRUE)

  obs_group <- h5_file[["obs"]]

  index_attr <- obs_group$attr_open("_index")
  index_col_name <- index_attr$read()
  index_attr$close()

  index <- obs_group[[index_col_name]][]

  all_obs_items <- obs_group$names
  data_cols <- setdiff(all_obs_items, c("__categories", index_col_name))

  obs_list <- list()

  if ("__categories" %in% all_obs_items) {
    categories_group <- obs_group[["__categories"]]
    category_names <- categories_group$names

    for (col in data_cols) {
      item <- obs_group[[col]]
      if (inherits(item, "H5D") && item$dims == length(index)) {
        if (col %in% category_names) {
          codes <- item[]
          levels <- categories_group[[col]][]
          codes_r <- codes + 1L
          codes_r[codes_r <= 0L] <- NA_integer_
          obs_list[[col]] <- factor(levels[codes_r], levels = levels)
        } else {
          obs_list[[col]] <- item[]
        }
      }
    }
  } else {
    for (col in data_cols) {
      item <- obs_group[[col]]
      if (inherits(item, "H5Group") &&
          all(c("codes", "categories") %in% item$names)) {
          codes <- item[["codes"]][]
          levels <- item[["categories"]][]
          codes_r <- codes + 1L
          codes_r[codes_r <= 0L] <- NA_integer_
          obs_list[[col]] <- factor(levels[codes_r], levels = levels)
      } else if (inherits(item, "H5D") && item$dims == length(index)) {
          obs_list[[col]] <- item[]
      }
    }
  }

  obs_df <- as.data.frame(obs_list, check.names = FALSE)
  rownames(obs_df) <- index

  return(obs_df)
}

remove_self_neighbors <- function(idx) {
  # Drop each cell's self-neighbor (value == row index); keep all if none.
  n <- nrow(idx); m <- ncol(idx)
  drop <- idx == seq_len(n)
  if (!any(drop)) return(idx)
  drop[rowSums(drop) == 0, m] <- TRUE
  matrix(t(idx)[t(!drop)], nrow = n, ncol = m - 1L, byrow = TRUE)
}

transfer_accuracy <- function(obs_df, orig_col, trans_col) {
  accuracy_df <- obs_df %>%
    mutate(
      orig = as.character(.data[[orig_col]]),
      trans = as.character(.data[[trans_col]])) %>%
    group_by(orig) %>%
    summarize(
      # na.rm matches Python (a null/NA prediction counts as incorrect but is
      # still in n_total) and prevents one NA from blanking the whole Total row
      n_correct = sum(orig == trans, na.rm = TRUE),
      n_total = n(),
      .groups = "drop") %>%
    bind_rows(
      summarize(
        .,
        orig = "Total",
        n_correct = sum(n_correct),
        n_total = sum(n_total))) %>%
    mutate(percent_correct = (n_correct / n_total) * 100) %>%
    rename(cell_type = orig)
  print(accuracy_df, n = Inf)
  return(accuracy_df)
}

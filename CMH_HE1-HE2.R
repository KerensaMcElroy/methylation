library(tidyverse)
library(data.table)
library(tidyr)
library(Deducer)
library(qqman)

# Try pooling by cohort.Sß
# Set working directory
setwd("/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_counts")

# Select chromosome (for testing)
# chr <- 29
# chr <- as.integer(chr)
# target_chrm_win <- "chrm29.1782210"


# Array on HE and a values
HE <- 1:2
a <- 1:14

# Functions used in for loop
countsToCases <- function(x, countcol = "count") {
  x %>%
    uncount(weights = !!sym(countcol), .remove = TRUE)
}

isEmpty <- function(x) {
  return(length(x)==0)   #insert 0 if CMH can't be computed
}

# main analysis function for each window
process_chrm_window <- function(target_win, merged_nonzero) {
  
  single_window_count <- merged_nonzero %>%
    filter(chrm_win == target_win) %>%
    pivot_longer(-chrm_win, names_to = c("a", "HE", "type"), names_sep = "_", values_to = "count") %>%
    mutate(count = as.integer(count))
  #remove stratum (e.g. ia) where counts are not available for both types (HE_1/HE_2 or NM/M counts) - less likely to encounter this if pooling counts across multiple animals
  
  filtered_window_count <- single_window_count %>%
    group_by(a, HE) %>%
    mutate(
      sum_type = sum(count[type == "meth" | type == "non"])
    ) %>%
    ungroup() %>%
    group_by(a, type) %>%
    mutate(
      sum_HE = sum(count[HE == "HE1" | HE == "HE2"])
    ) %>%
    ungroup() %>%
    group_by(a) %>%
    filter(all(sum_type != 0) && all(sum_HE != 0)) %>%
    ungroup() %>%
    dplyr::select(-sum_type, -sum_HE)
  
  filtered_window_cases <- countsToCases(filtered_window_count)

  num_strata <- filtered_window_cases %>%
    distinct(a) %>%
    n_distinct()
  
  if (num_strata < 2) {
    cat("Only one stratum is available. Skipping CMH test.\n")
    return(data.frame(chrm_win = target_win, CMH_st = 0, CMH_p = 0, CMH_estimate = 0))
  }
  
  tryCatch({
    test <- contingency.tables(
      row.vars = type,
      col.vars = HE,
      stratum.var = a,
      data = filtered_window_cases
    )
    print(test)
    test <- add.mantel.haenszel(test)
  }, 
  warning = function(w) {
    # Handling the warning
    # Print the warning message and chrm_win value
    warning(paste("Warning: test failed in chrm_win:", chrm_win))
    warning(w)
  })
  
  CMH_st <- as.numeric(test$"type by HE"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$statistic)
  if (isEmpty(CMH_st)){
    CMH_st = 0
  }
  
  CMH_p <- (test$"type by HE"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$p.value)
  if (isEmpty(CMH_p)){
    CMH_p = 0
  }
  
  CMH_estimate <- as.numeric(test$"type by HE"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$estimate)
  if (isEmpty(CMH_estimate)){
    CMH_estimate = 0
  }
  return(data.frame(chrm_win = target_win, CMH_st = CMH_st, CMH_p = CMH_p, CMH_estimate = CMH_estimate))
}

#process each chromosome in turn
for (chr in 1:29) {
  # Read files and merge data
  merged_data <- expand_grid(a = a, HE=HE) %>%
    mutate(file_name = paste0("window_sum.count_HE", HE, "_chrm", chr, "_ia", a, ".txt")) %>%
    rowwise() %>%
    mutate(data = list(fread(file.path(getwd(), file_name), header = FALSE))) %>%
    ungroup() %>%
    mutate(data = map(data, ~ setNames(.x, c("non", "total")))) %>%
    unnest(data) %>%
    group_by(file_name) %>%
    mutate(row_number = row_number()) %>%
    ungroup() %>%
    mutate(meth = total - non) %>%
    dplyr::select(-total, -file_name) %>%
    pivot_longer(cols = c(non, meth), names_to = "type", values_to = "value") %>%
    pivot_wider(names_from=c(a, HE, type), values_from = value) %>%
    mutate(chrm_win = paste0("chrm", chr, ".", row_number)) %>%
    dplyr::select(-row_number) %>%
    rename_with(~paste0("a", str_extract(.x, "\\d+"), "_", "HE", str_extract(.x, "(?<=_)\\d"), "_", str_extract(.x, "(?<=_)\\D+$")), -chrm_win) %>%
    dplyr::select(chrm_win, everything())
  
  #remove rows with only zero entries
  merged_nonzero <- merged_data %>%
    filter(rowSums(across(-chrm_win, ~ . == 0)) != (ncol(.)-1)) 
  
  results <- map_df(merged_nonzero$chrm_win, process_chrm_window, merged_nonzero = merged_nonzero)
  
  write.csv(results, paste0("~/R/methylation/results/tables/cmh_", chr,".csv"), row.names = FALSE)
  
  filtered_results <- results %>%
    mutate(SNP=chrm_win) %>%
    separate(chrm_win, into = c("chrm", "bp"), sep = "\\.", convert = TRUE) %>%
    mutate(chrm = as.numeric(gsub("[^0-9]", "", chrm))) %>%
    filter(CMH_st > 0) %>%
    mutate(corrected = p.adjust(CMH_p, method = "bonferroni")) %>%
    mutate(CMH_estimate = na_if(CMH_estimate, Inf))
  
  # Find significant windows
  highlighted_windows <- filtered_results %>%
    filter(corrected < 0.05) %>%
    pull(SNP)
  
  # manhattan plots for CMH statistic, estimate, and p-values
  pdf(file = paste0("~/R/methylation/results/figures/man_st_", chr, ".pdf"), width = 8, height = 5)
  manhattan_st <- manhattan(
    filtered_results,
    chr = "chrm",  # Column containing chromosome information
    bp = "bp",  # Column containing base pair positions (if available, or set to NULL)
    p = "CMH_st",  # Column containing CMH_p values
    logp=FALSE,
    ylim = c(0, max(filtered_results$CMH_st)*1.1),
    suggestiveline = FALSE,
    genomewideline = FALSE,
    highlight = highlighted_windows,
    main = "CMH_st",
    ylab = "CMH statistic"
  )
  dev.off()
  
  pdf(file = paste0("~/R/methylation/results/figures/man_est_", chr, ".pdf"), width = 8, height = 5)
  manhattan_est <- manhattan(
    filtered_results,
    chr = "chrm",  # Column containing chromosome information
    bp = "bp",  # Column containing base pair positions (if available, or set to NULL)
    p = "CMH_estimate",  # Column containing CMH_p values
    logp=FALSE,
    ylim = c(0, max(filtered_results$CMH_estimate, na.rm=TRUE)*1.1),
    suggestiveline = FALSE,
    genomewideline = FALSE,
    highlight = highlighted_windows,
    main = "CMH_estimate",
    ylab = "odds ratio estimate"
  )
  dev.off()
  
  pdf(file = paste0("~/R/methylation/results/figures/man_p_", chr, ".pdf"), width = 8, height = 5)
  manhattan_p <- manhattan(
    filtered_results,
    chr = "chrm",  # Column containing chromosome information
    bp = "bp",  # Column containing base pair positions (if available, or set to NULL)
    p = "CMH_p",  # Column containing CMH_p values
    ylim = c(0, 30),
    suggestiveline = FALSE,
    genomewideline = FALSE,
    highlight = highlighted_windows,
    main = "raw p-values"
  )
  dev.off()
}





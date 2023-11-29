library(tidyverse)
library(data.table)
library(tidyr)
library(Deducer)
library(qqman)
library(readxl)
library(janitor)


# Try pooling by cohort.Sß
# Set working directory
setwd("/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_counts")

# Select chromosome (for testing)
# chr <- 29
# chr <- as.integer(chr)
target_chrm_win <- "chrm_28.1782210"


# Array on HE and a values
HE <- 1:2
a <- 1:14
cohort <- 7:8
room <- 1:4
pool <- "cohort"


# Functions used in for loop
countsToCases <- function(x, countcol = "count") {
  x %>%
    uncount(weights = !!sym(countcol), .remove = TRUE) %>%
  print()
}

isEmpty <- function(x) {
  return(length(x)==0)   #insert 0 if CMH can't be computed
}

# additional analysis to set room or cohort as level in CMH analysis

key <- read_csv("~/R/methylation/data/sample_key.csv")
data_path <- "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/biochem/CC7+8  biochem for Yutao 05_02_21.xlsx"
features <- data_path %>%
  read_excel(col_names = FALSE) %>% 
  dplyr::select(-1:-2) %>% 
  slice(1:3) %>%
  t() %>%
  as_tibble() %>%
  rename('ID' = 1, 'cohort' = 2, 'room' = 3) %>%
  mutate_if(is.numeric, as.character)

animal_map <- key %>%
  dplyr::select(`window mean ID`, `Animal ID`) %>%
  distinct() %>%
  drop_na() %>%
  clean_names() %>%
  mutate(window_mean_id = as.character(window_mean_id),
         animal_id = as.character(animal_id)) %>%
  full_join(features, join_by(animal_id == ID)) 



# main analysis function for each window
 
process_chrm_window <- function(target_win, merged_nonzero) {
  single_window <- merged_nonzero %>%
    filter(chrm_win == target_win)
  single_window_count <- single_window %>%
    pivot_longer(c(-chrm_win, -chr, -winMin, -winMax), names_to = c("pool", "HE", "type"), names_sep = "_", values_to = "count") %>%
    mutate(count = as.integer(count))
  #remove stratum (e.g. ia) where counts are not available for both types (HE_1/HE_2 or NM/M counts) - less likely to encounter this if pooling counts across multiple animals
  
  filtered_window_count <- single_window_count %>%
    group_by(pool, HE) %>%
    mutate(
      sum_type = sum(count[type == "meth" | type == "non"])
    ) %>%
    ungroup() %>%
    group_by(pool, type) %>%
    mutate(
      sum_HE = sum(count[HE == "HE1" | HE == "HE2"])
    ) %>%
    ungroup() %>%
    group_by(pool) %>%
    filter(all(sum_type != 0) && all(sum_HE != 0)) %>%
    ungroup() %>%
    dplyr::select(-sum_type, -sum_HE) %>%
    print()

  filtered_window_cases <- countsToCases(filtered_window_count)
  num_strata <- filtered_window_cases %>%
    distinct(pool) %>%
    n_distinct() 

  if (num_strata < 2) {
    return(data.frame(chrm_win = target_win, CMH_st = 0, CMH_p = 0, CMH_estimate = 0, counts = 0, winMin= single_window$winMin, winMax= single_window$winMax, chr= single_window$chr))
  }
  
  tryCatch({
    test <- contingency.tables(
      row.vars = type,
      col.vars = HE,
      stratum.var = pool,
      data = filtered_window_cases
    )
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
  
 # print(sapply(test$"type by HE", function(a_value) a_value$total))
  total_mean <- mean(sapply(test$"type by HE", function(pool_value) if("total" %in% names(pool_value) && is.numeric(pool_value$total)) pool_value$total else 0))
  
  # Print the result
  return(data.frame(chrm_win = target_win, CMH_st = CMH_st, CMH_p = CMH_p, CMH_estimate = CMH_estimate, counts = total_mean, winMin= single_window$winMin, winMax= single_window$winMax, chr= single_window$chr))
}

#process each chromosome in turn
for (chrm in 28:28) {
  # Read files and merge data
  win_pos <- read_csv(paste0("/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_means/window_positions_chrm", chrm, ".csv"))
  merged_data <- expand_grid(a = a, HE=HE) %>%
    mutate(file_name = paste0("window_sum.count_HE", HE, "_chrm", chrm, "_ia", a, ".txt")) %>%
    rowwise() %>%
    mutate(data = list(fread(file.path(getwd(), file_name), header = FALSE))) %>%
    ungroup() %>%
    mutate(data = map(data, ~ setNames(.x, c("non", "total")))) %>%
    unnest(data) %>%
    group_by(file_name) %>%
    mutate(row_number = row_number()) %>%
    ungroup() %>%
    mutate(meth = total - non) %>%
    mutate(animal = as.character(a)) %>%
    left_join(animal_map, join_by(animal == window_mean_id)) %>%
    dplyr::select(-total, -file_name) %>%
    group_by(!!sym(pool), HE, row_number) %>%
    summarise(across(c(non, meth), \(x) sum(x, na.rm = TRUE))) %>%
    pivot_longer(cols = c(non, meth), names_to = "type", values_to = "value") %>%
    pivot_wider(names_from=c(!!sym(pool), HE, type), values_from = value) %>%
    mutate(chrm_win = paste0("chrm_", chrm, ".", row_number)) %>%
    left_join(win_pos, join_by(chrm_win == chrm.window)) %>%
    mutate(chr = as.character(chrm)) %>%
    dplyr::select(-row_number) %>%
    rename_with(~paste0(pool, str_extract(.x, "\\d+"), "_", "HE", str_extract(.x, "(?<=_)\\d"), "_", str_extract(.x, "(?<=_)\\D+$")), c(-chrm_win, -winMin, -winMax, -chr)) %>%
    dplyr::select(chrm_win, chr, winMin, winMax, everything()) 

  #remove rows with only zero entries
  merged_nonzero <- merged_data %>%
    filter(rowSums(across(c(-chrm_win, -winMin, -winMax, -chr), ~ . == 0)) != (ncol(.)-1)) 
  
  results <- map_df(merged_nonzero[0:1000,]$chrm_win, process_chrm_window, merged_nonzero = merged_nonzero)
  
  write.csv(results, paste0("~/R/methylation/results/tables/cmh_",  chrm, "_", pool, ".csv"), row.names = FALSE)
  
  filtered_results <- results %>%
    mutate(SNP=chrm_win) %>%
    separate(chrm_win, into = c("chrm", "index"), sep = "\\.", convert = TRUE) %>%
    mutate(chrm = as.numeric(gsub("[^0-9]", "", chrm))) %>%
    filter(CMH_st > 0) %>%
    mutate(corrected = p.adjust(CMH_p, method = "bonferroni")) %>%
    mutate(CMH_estimate = na_if(CMH_estimate, Inf))
  
  # Find significant windows
  highlighted_windows <- filtered_results %>%
    filter(corrected < 0.05) %>%
    pull(SNP)
  
  # manhattan plots for CMH statistic, estimate, and p-values
  pdf(file = paste0("~/R/methylation/results/figures/", chrm, "_", pool, "_man_st.pdf"), width = 8, height = 5)
  manhattan_st <- manhattan(
    filtered_results,
    chr = "chrm",  # Column containing chromosome information
    bp = "winMin",  # Column containing base pair positions (if available, or set to NULL)
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
  
  pdf(file = paste0("~/R/methylation/results/figures/", chrm, "_", pool, "_man_est.pdf"), width = 8, height = 5)
  manhattan_est <- manhattan(
    filtered_results,
    chr = "chrm",  # Column containing chromosome information
    bp = "winMin",  # Column containing base pair positions (if available, or set to NULL)
    p = "CMH_estimate",  # Column containing CMH_p values
    logp=TRUE,
    ylim = c(-2, 2),
    suggestiveline = 0,
    genomewideline = FALSE,
    highlight = highlighted_windows,
    main = "CMH_estimate",
    ylab = "log odds ratio estimate"
  )
  dev.off()
  
  pdf(file = paste0("~/R/methylation/results/figures/", chrm, "_", pool, "_man_p.pdf"), width = 8, height = 5)
  manhattan_p <- manhattan(
    filtered_results,
    chr = "chrm",  # Column containing chromosome information
    bp = "winMin",  # Column containing base pair positions (if available, or set to NULL)
    p = "CMH_p",  # Column containing CMH_p values
    ylim = c(0, 30),
    suggestiveline = FALSE,
    genomewideline = FALSE,
    highlight = highlighted_windows,
    main = "raw p-values"
  )
  dev.off()
  
  scat_bp <- ggplot(filtered_results, aes(x = winMin, y = counts, colour = corrected < 0.05)) +
    geom_point(aes(alpha = !corrected >= 0.05), size = 2) +
    scale_colour_manual(values = c("TRUE" = "green", "FALSE" = "black"), guide = FALSE) +
    guides(alpha = FALSE) 
  ggsave(filename= paste0("~/R/methylation/results/figures/", chrm, "_", pool, "_scatter_bp.pdf"),plot=scat_bp, width = 10, height = 15)
  scat_st <- ggplot(filtered_results, aes(x = CMH_st, y = counts, colour = corrected < 0.05)) +
    geom_point(aes(alpha = !corrected >= 0.05), size = 2) +
    scale_colour_manual(values = c("TRUE" = "green", "FALSE" = "black")) +
    guides(alpha = FALSE) 
  ggsave(filename= paste0("~/R/methylation/results/figures/", chrm, "_", pool, "_scatter_st.pdf"),plot=scat_st, width = 10, height = 10)
}

#todo - move plots to separate script
#set up permutation
#set up whole genome
#set up sbatch job.


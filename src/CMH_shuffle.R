library(tidyverse)
library(data.table)
library(tidyr)
library(Deducer)
library(qqman)
library(readxl)
library(janitor)

# Set working directory
setwd("/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_counts")

# Options for testing
# chr <- as.integer(28)
# target_chrm_win <- "chrm_28.1782210"


# Array on HE and a values
HE <- 1:2
a <- 1:14
cohort <- 7:8
room <- 1:4
pool <- "animal" #one of "cohort", "room", "animal". Determins the stratification for CMH test.

# add room and cohort metadata to methylation data. 

# contains room and cohort info
data_path <- "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/biochem/CC7+8  biochem for Yutao 05_02_21.xlsx"

metadata <- data_path %>%
  read_excel(col_names = FALSE) %>% 
  dplyr::select(-1:-2) %>% 
  slice(1:3) %>%
  t() %>%
  as_tibble() %>%
  rename('id' = 1, 'cohort' = 2, 'room' = 3) %>%
  mutate_if(is.numeric, as.character)

# contains mappings between animal IDs in methylation and biochem data.
key <- read_csv("~/R/methylation/data/sample_key.csv")

animal_map <- key %>%
  dplyr::select(`window mean ID`, `Animal ID`) %>%
  distinct() %>%
  drop_na() %>%
  clean_names() %>%
  rename(id=animal_id)%>%
  mutate(window_mean_id = as.character(window_mean_id),
         id = as.character(id))

metadata <- animal_map %>%
  full_join(metadata, join_by(id == id)) 

# Helper function to filter and process counts
window_counts <- function(single_window) {
  single_window %>%
    pivot_longer(c(-chrm_win, -chr, -winMin, -winMax), 
                 names_to = c("pool", "HE", "type"), 
                 names_sep = "_", 
                 values_to = "count") %>%
    mutate(count = as.integer(count))
}

# Helper function thats filters counts to remove windows incompatible with CMH due to zeros.
filter_counts <- function(win_counts) {
  win_counts %>%
    group_by(pool, HE) %>%
    mutate(sum_type = sum(count[type %in% c("meth", "non")])) %>%
    ungroup() %>%
    group_by(pool, type) %>%
    mutate(sum_HE = sum(count[HE %in% c("HE1", "HE2")])) %>%
    ungroup() %>%
    group_by(pool) %>%
    filter(all(sum_type != 0) && all(sum_HE != 0)) %>%
    ungroup() %>%
    dplyr::select(-sum_type, -sum_HE)
}

# Helper function to convert counts to cases format
counts_to_cases <- function(x, countcol = "count") {
  x %>%
    uncount(weights = !!sym(countcol), .remove = TRUE) 
}

# main analysis function for each window
process_chrm_window <- function(target_win, merged_nonzero) {
  single_window <- merged_nonzero %>%
    filter(chrm_win == target_win)
  
  filtered_window_cases <- window_counts(single_window) %>%
    filter_counts() %>%
    counts_to_cases()
  
  total_mean <- mean(sapply(test$"type by HE", function(pool_value) if("total" %in% names(pool_value) && is.numeric(pool_value$total)) pool_value$total else 0))
  
  # Define default values
  default_values <- c(
    chrm_win = target_win,
    CMH_st = 0,
    CMH_p = 0,
    CMH_estimate = 0,
    mean_counts = total_mean,
    winMin = single_window$winMin,
    winMax = single_window$winMax,
    chr = single_window$chr
  )
  
  # check we have at least two strata and return without calculating CMH otherwise
  if (num_strata < 2) {
    return(as.data.frame(default_values))
  }
  
  #run the CMH test
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
    warning(paste("Warning: test failed in chrm_win:", target_win))
    warning(w)
  })
  
  #extract statistics from data structure
  CMH_st <- coalesce(as.numeric(test$"type by HE"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$statistic), 0)
  CMH_p <- coalesce(test$"type by HE"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$p.value, 0)
  CMH_estimate <- coalesce(as.numeric(test$"type by HE"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$estimate), 0)
  
  return(as.data.frame(c(default_values, CMH_st = CMH_st, CMH_p = CMH_p, CMH_estimate = CMH_estimate)))
}

# Define a function to read and merge files for a single chromosome
merge_chrm <- function(a, HE, chrm, metadata) {
  path <- "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_means/window_positions_chrm"
  win_pos <- read_csv(paste0(path, chrm, ".csv"))
  merged_data <- expand_grid(a = a, HE=HE) %>%
    mutate(file_name = paste0("window_sum.count_HE", HE, "_chrm", chrm, "_ia", a, ".txt")) %>%
    rowwise() %>%
    mutate(data = list(fread(file.path(getwd(), file_name), header = FALSE))) %>%
    ungroup() %>%
    mutate(data = map(data, ~ setNames(.x, c("non", "total")))) %>%
    unnest(data) %>%
    group_by(file_name) %>%
    mutate(win_num = row_number()) %>%
    ungroup() %>%
    mutate(meth = total - non) %>%
    mutate(win_animal = as.character(a)) %>%
    mutate(HE = as.character(HE)) %>%
    mutate(win_num = as.character(win_num)) %>%
    dplyr::select(win_animal, HE, win_num, non, meth) %>%
    left_join(metadata, join_by(win_animal == window_mean_id)) 
  return(merged_data)
}


shuffle_data <- function(merge_data) {
  shuffle <- sample(nrow(merge_data)) 
  shuffled_data <- chrm_data %>%
    mutate(
      non = non[shuffle],
      meth = meth[shuffle]
    )
}

#process each chromosome in turn
for (chrm in 1:29) {
  # Read files and merge data

    group_by(!!sym(pool), HE, row_number) %>%
    summarise(across(c(non, meth), \(x) sum(x, na.rm = TRUE))) %>%
    print() %>%
    pivot_longer(cols = c(non, meth), names_to = "type", values_to = "value") %>%
    pivot_wider(names_from=c(!!sym(pool), HE, type), values_from = value) %>%
    mutate(chrm_win = paste0("chrm_", chrm, ".", row_number)) %>%
    left_join(win_pos, join_by(chrm_win == chrm.window)) %>%
    mutate(chr = as.character(chrm)) %>%
    dplyr::select(-row_number) %>%
    print() %>%
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


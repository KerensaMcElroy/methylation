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

 chr <- as.integer(28)
 target_chrm_win <- "chrm_28.1782210"


# Array on HE and a values
HE <- 1:2
a <- 1:14
cohort <- 7:8
room <- 1:4
pool <- "id" #one of "cohort", "room", "id". Determines the stratification for CMH test.

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
  left_join(metadata, join_by(id == id)) 

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

  # Define default values
  default_values <- data.frame(
    chrm_win = target_win,
    CMH_st = 0,
    CMH_p = 0,
    CMH_estimate = 0,
    mean_counts = 0,
    winMin = single_window$winMin,
    winMax = single_window$winMax,
    chr = single_window$chr
  )
  # check we have at least two strata and return without calculating CMH otherwise
  num_strata <- length(unique(filtered_window_cases$pool))
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
  total_mean <- mean(sapply(test$"type by HE", function(pool_value) if("total" %in% names(pool_value) && is.numeric(pool_value$total)) pool_value$total else 0),)
  
  results <- default_values %>%
    mutate(CMH_st = CMH_st, CMH_p = CMH_p, CMH_estimate = CMH_estimate, mean_counts = total_mean)
  return(results)
}

# Define a function to read and merge files for a single chromosome
merge_chrm <- function(a, HE, chrm, metadata) {
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

#helper function to randomly shuffle data. Keeps meth and non-meth together and shuffles relative to everything else.
shuffle_data <- function(merge_data) {
  shuffle <- sample(nrow(merge_data)) 
  shuffled_data <- chrm_data %>%
    mutate(
      non = non[shuffle],
      meth = meth[shuffle]
    )
  return(shuffled_data)
}

#helper function to pool based on stratification (id, room, cohort)
#to do ... move some of the metadata joining out of here.

pool_data <- function(data) {
  path <- "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_means/window_positions_chrm"
  win_pos <- read_csv(paste0(path, chrm, ".csv"))
  pooled <- data %>%
    group_by(!!sym(pool), HE, win_num) %>%
    summarise(across(c(non, meth), \(x) sum(x, na.rm = TRUE))) %>%
    pivot_longer(cols = c(non, meth), names_to = "type", values_to = "value") %>%
    pivot_wider(names_from=c(!!sym(pool), HE, type), values_from = value) %>%
    mutate(chrm_win = paste0("chrm_", chrm, ".", win_num)) %>%
    left_join(win_pos, join_by(chrm_win == chrm.window)) %>%
    mutate(chr = as.character(chrm)) %>%
    dplyr::select(-win_num) %>%
    rename_with(~paste0(pool, str_extract(.x, "\\d+"), "_", "HE", str_extract(.x, "(?<=_)\\d"), "_", str_extract(.x, "(?<=_)\\D+$")), c(-chrm_win, -winMin, -winMax, -chr)) %>%
    dplyr::select(chrm_win, chr, winMin, winMax, everything()) 
  return(pooled)
}
  

#process each chromosome in turn
for (chrm in 1:9) {
  # Read files and merge data
  chrm_data <- merge_chrm(a, HE, chrm, metadata) 
  chrm_shuf <- shuffle_data(chrm_data) 
  chrm_pool <- pool_data(chrm_shuf) 

  #remove rows with only zero entries. Needs to be done here due to pooling.
  columns_to_drop <- c("chrm_win", "winMin", "winMax", "chr")
  chrm_nonzero <- chrm_pool %>%
    filter(rowSums(across(-columns_to_drop, ~ . != 0)) > 0) 

  results <- map_df(chrm_nonzero$chrm_win[1:100], process_chrm_window, merged_nonzero = chrm_nonzero) %>%
    arrange(winMin)
  
  
  write.csv(results, paste0("~/R/methylation/results/tables/cmh_",  chrm, "_", pool, ".csv"), row.names = FALSE)
}

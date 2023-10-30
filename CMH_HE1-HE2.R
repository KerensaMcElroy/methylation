library(tidyverse)
library(data.table)
library(tidyr)
library(Deducer)

# Try pooling by cohort.Sß
# Set working directory
setwd("/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_counts")

# Select chromosome
chr <- 29
chr <- as.integer(chr)
target_chrm_win <- "chrm29.390"


# Array on HE and a values
HE <- 1:2
a <- 1:14


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
  
# question for shannon ... is the chrm window just a sequential number or the start position? Can't work it out from your code.

# Print the merged data
print(merged_data)

#remove rows with only zero entries
merged_nonzero <- merged_data %>%
  filter(rowSums(across(-chrm_win, ~ . == 0)) != (ncol(.)-1)) 

single_window_count <- merged_nonzero %>%
  filter(chrm_win == target_chrm_win) %>%
  pivot_longer(-chrm_win, names_to = c("a", "HE", "type"), names_sep = "_", values_to = "count") %>%
  mutate(count = as.integer(count))

# to do - summarise by combination of room and cohort. Will need to process from biochem and sample_key. Also try summarising only by cohort.
# doesn't make sense to group like this by a1-a7 and a8-a14
summarized_window_count <- single_window_count %>%
  mutate(a_group = ifelse(as.integer(substring(a, 2)) <= 7, "a1-a7", "a8-a14")) %>%
  group_by(HE, type, a_group) %>%
  summarise(sum_count = sum(count), .groups = "drop") %>%
  group_by(a_group) %>%
  filter(!sum(sum_count == 0)) %>%
  ungroup()

#remove stratum (e.g. ia) where counts are not available for both types (HE_1/HE_2 or NM/M counts) - less likely to encounter this if pooling counts across multiple animals
# Filter rows based on the conditions

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

 
#check with shannon that I've understood this filtering correctly
#question for shannon - should we do this, or should we set the zeros to one? In effect we're throwing away any windows where the chi-squared goes to infinity, and therefore potentially the most significant results.

#grouping may be problematic due to batch effect

countsToCases<- function(x, countcol = "count") {
  x %>%
    uncount(weights = !!sym(countcol), .remove = TRUE)
}

filtered_window_cases <- countsToCases(filtered_window_count)
summarized_window_cases <- countsToCases(summarized_window_count, countcol = "sum_count")

num_strata <- summarized_window_cases %>%
  distinct(a_group) %>%
  n_distinct()

tryCatch({
  test <- contingency.tables(
    row.vars = type,
    col.vars = HE,
    stratum.var = a_group,
    data = summarized_window_cases
  )
  
  test <- add.chi.squared(test)
  test <- add.likelihood.ratio(test)
  test <- add.mantel.haenszel(test)
  
}, warning = function(w) {
  # Handling the warning
  # Print the warning message and chrm_win value
  warning(paste("Warning: test failed in chrm_win:", chrm_win))
  warning(w)
})
  
  if (inherits(test, "try-error")) {
    # Error occurred, assign 0 to chi_squared and p_value for all strata
    chi_squared <- rep(0, num_strata)
    p_value <- rep(0, num_strata)
  } else {
    # Extract chi-squared and p-value for each stratum
    chi_squared <- vector("numeric", length = num_strata)
    p_value <- vector("numeric", length = num_strata)
    
    for (i in 1:num_strata) {
      tryCatch({
        chi_squared[i] <- as.numeric(test$`type by ENV`$tests$`Chi Squared`[i][[1]]$asymptotic$statistic)
        p_value[i] <- test$`type by ENV`$tests$`Chi Squared`[i][[1]]$asymptotic$p.value
        
      }, error = function(e) {
        # Handling the error for individual stratum
        # Assign 0 to chi_squared and p_value for the current stratum
        chi_squared[i] <- 0
        p_value[i] <- 0
      })
    }
  }
  
  result <- list(Chisq_ia = chi_squared, Chisq_p_ia = p_value)

  
#insert 0 if CMH can't be computed
isEmpty <- function(x) {
   return(length(x)==0)
}
  
CMH_st <- as.numeric(test$"type by ENV"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$statistic)
if (isEmpty(CMH_st)){
  CMH_st = 0
}
CMH_p <- (test$"type by ENV"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$p.value)
if (isEmpty(CMH_p)){
  CMH_p = 0
}
CMH_estimate <- as.numeric(test$"type by ENV"$cross.strata.tests$"Mantel-Haenszel"[1][[1]]$asymptotic$estimate)
if (isEmpty(CMH_estimate)){
  CMH_estimate = 0
}
  
# View the filtered tibble
#results = data.frame(matrix(nrow = nrow(merged_nonzero),ncol=32), stringsAsFactors=FALSE)
#names(results) = c("chrm.window", " Chisq_ia1", " Chisq_p_ia1", " Chisq_ia2", " Chisq_p_ia2", " Chisq_ia3", " Chisq_p_ia3", " Chisq_ia4", " Chisq_p_ia4", " Chisq_ia5", " Chisq_p_ia5", " Chisq_ia6", " Chisq_p_ia6", " Chisq_ia7", " Chisq_p_ia7", " Chisq_ia8", " Chisq_p_ia8", " Chisq_ia9", " Chisq_p_ia9", " Chisq_ia10", " Chisq_p_ia10", " Chisq_ia11", " Chisq_p_ia11", " Chisq_ia12", " Chisq_p_ia12", " Chisq_ia13", " Chisq_p_ia13", " Chisq_ia14", " Chisq_p_ia14", " CMH_st", " CMH_p", " CMH_estimate")

results = data.frame(matrix(nrow = dimJoin[1],ncol=8), stringsAsFactors=FALSE)
names(results) = c("chrm.window", " Chisq_ia1", " Chisq_p_ia1", " Chisq_ia2", " Chisq_p_ia2", " CMH_st", " CMH_p", " CMH_estimate")


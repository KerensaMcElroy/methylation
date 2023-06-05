library(tidyverse)
library(data.table)
library(tidyr)
library(Deducer)

# Set working directory
setwd("/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_counts")

# Select chromosome
chr <- 29
chr <- as.integer(chr)

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
  filter(chrm_win == "chrm29.2") %>%
  pivot_longer(-chrm_win, names_to = c("a", "HE", "type"), names_sep = "_", values_to = "count") %>%
  mutate(count = as.integer(count))

summarized_window_count <- single_window_count %>%
  mutate(a_group = ifelse(as.integer(substring(a, 2)) <= 7, "a1-a7", "a8-a14")) %>%
  group_by(HE, type, a_group) %>%
  summarise(sum_count = sum(count))

filtered_window_count_pool <- single_window_count_pool(merged_nonzero)


#remove stratum (e.g. ia) where counts are not available for both types (HE_1/HE_2 or NM/M counts) - less likely to encounter this if pooling counts across multiple animals
# Filter rows based on the conditions
filtered_window_count <- single_window_count %>%
  group_by(a) %>%
  filter(!any(count == 0)) %>%
  ungroup()
#check with shannon that I've understood this filtering correctly

#shannons next step is to remove the counts? (counts to cases)
# get number of remaining animals...

countsToCases<- function(x, countcol = "count") {
  x %>%
    uncount(weights = !!sym(countcol), .remove = TRUE)
}

filtered_window_cases <- countsToCases(filtered_window_count)
summarized_window_cases <- countsToCases(summarized_window_count, countcol = "sum_count")
test <-contingency.tables(
  row.vars=type,
  col.vars=HE,
  stratum.var=a_group,data=summarized_window_cases)
test <- add.chi.squared(test)
test <- add.likelihood.ratio(test)
test <- add.mantel.haenszel(test)

# View the filtered tibble
results = data.frame(matrix(nrow = nrow(merged_nonzero),ncol=32), stringsAsFactors=FALSE)
names(results) = c("chrm.window", " Chisq_ia1", " Chisq_p_ia1", " Chisq_ia2", " Chisq_p_ia2", " Chisq_ia3", " Chisq_p_ia3", " Chisq_ia4", " Chisq_p_ia4", " Chisq_ia5", " Chisq_p_ia5", " Chisq_ia6", " Chisq_p_ia6", " Chisq_ia7", " Chisq_p_ia7", " Chisq_ia8", " Chisq_p_ia8", " Chisq_ia9", " Chisq_p_ia9", " Chisq_ia10", " Chisq_p_ia10", " Chisq_ia11", " Chisq_p_ia11", " Chisq_ia12", " Chisq_p_ia12", " Chisq_ia13", " Chisq_p_ia13", " Chisq_ia14", " Chisq_p_ia14", " CMH_st", " CMH_p", " CMH_estimate")

perform_contingency_tests <- function(counts_matrix, temp_data_transposed) {
  test <- tryCatch({
    contingency_tables(
      row.vars = d(type),
      col.vars = d(ENV),
      stratum.var = ia,
      data = temp_data_transposed
    ) %>%
      add.chi.squared() %>%
      add.likelihood.ratio() %>%
      add.mantel.haenszel()
  }, error = function(e) {
    NULL
  })
  
  if (is.null(test)) {
    return(list(Chisq_ia = 0, Chisq_p_ia = 0))
  }
  
  p_values <- map_dbl(test$"type by ENV"$tests$"Chi Squared", ~ .x$asymptotic$p.value)
  
  list(
    Chisq_ia = map_dbl(test$"type by ENV"$tests$"Chi Squared", ~ .x$asymptotic$statistic),
    Chisq_p_ia = ifelse(is.nan(p_values), 0, p_values)
  )
}

# Convert counts to cases
countsToCases <- function(x, countcol = "counts") {
  idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
  x[[countcol]] <- NULL
  x[idx, ]
}

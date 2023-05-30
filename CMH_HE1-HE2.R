library(tidyverse)
library(data.table)

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

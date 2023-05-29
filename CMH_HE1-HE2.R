library(tidyverse)
library(data.table)

# Set working directory
setwd("/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_counts")

# Select chromosome
chr <- 29
chr <- as.integer(chr)

# Array on HE and a values
HE <- 1:3
a <- 1:14

# Read files and merge data
merged_data <- expand_grid(HE = HE, a = a) %>%
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
  pivot_wider(names_from=c(HE, a, type), values_from = value) %>%
  dplyr::select(-row_number) 




# Print the merged data
print(merged_data)

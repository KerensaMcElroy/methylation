### Load libraries
library(conflicted)
library(tidyverse)
library(readxl)
library(mixOmics)
library(janitor)

# Path to data
data_path <- "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/biochem/CC7+8  biochem for Yutao 05_02_21.xlsx"
#We need to turn all this into a tidy data frame. Columns need to be ID, room, cohort, day, treatment, and biochem features. Starting by using the cheat sheet to get ID, cohort, and room. 

features <- data_path %>%
  read_excel(col_names = FALSE) %>% 
  dplyr::select(-1:-2) %>% 
  slice(1:3) %>%
  t() %>%
  as_tibble() %>%
  rename('ID' = 1, 'cohort' = 2, 'room' = 3) %>%
  mutate_if(is.numeric, as.character)
print(features)

#Now, we read in the first row of a single biochem sheet so that we can get the column names. This is necessary because due to the formatting of the excel sheets, we need to exclude the first few rows, which discards the column names.
col_names <- read_excel(path=data_path, sheet=2,range="A1:Y1",col_names=TRUE)
#Now we: * select the data range on each sheet and use map_df to combine into a single data frame, discard 'cheat sheet' from the list of names, use mutate to add a treatment column, use pivot and join functions to format resulting tibble
biochem <- data_path %>%
  excel_sheets() %>%
  set_names() %>% 
  keep(names(.)!="cheat sheet") %>%
  map_df(~ read_excel(path = data_path, sheet = .x, range =  anchored("A4", c(13, 25)), na='.', col_names=names(col_names)), .id = "biochem") %>%
  rename('days' = 'ID') %>%
  mutate(treatment = case_when(
    days > 17 ~ "pens",
    days > 12 ~ "recovery",
    days > 5  ~ "hot",
    days <=5 ~ "pre"
  )) %>%
  pivot_longer(cols = 3:26, names_to = 'ID') %>%
  full_join(features, by = 'ID') %>%
  pivot_wider(names_from = 'biochem', values_from = 'value')
print(biochem)

biochem_meth <- biochem %>%
  dplyr::filter(days %in% c(5, 12, 38))

# connect to methylation database
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_means/meth_db.sqlite")
# select data for given chromosome
query <- "SELECT methylation, ID, treatment,winMin FROM your_table WHERE chromosome = 29"

pca_data <- dbGetQuery(con, query) %>%
  as_tibble()

dbDisconnect(con)

pca_wide <- pca_data %>% pivot_wider(values_from = "methylation", names_from = "winMin")
pca_wide <- pca_wide %>% 
  dplyr::select(where(~mean(is.na(.)) < .2))
var_df <- pca_wide %>% 
  dplyr::select(-ID,-treatment) %>%
  summarise(across(.cols = everything(), ~var(.x, na.rm = TRUE))) %>% 
  pivot_longer(cols = everything()) %>% 
  arrange(desc(value)) %>%
  mutate(order = 1:n())

var_1000 = var_df %>%
  slice_max(n = 1000, order_by = value)

hist_plot <- var_df %>% 
  mutate(included= case_when(order <= 1000 ~ TRUE, .default = FALSE)) %>%
  ggplot(aes(x = value, fill = included)) +
  geom_histogram(bins=100)  +
  scale_fill_brewer(palette="Accent")

hist_plot

for_pca <- pca_wide %>% 
  dplyr::select(var_1000$name)
tidy_pca <- as_tibble(result.pca.multi[["x"]]) %>%
  add_column(ID = pca_wide$ID, treatment = pca_wide$treatment) %>%
  inner_join(biochem_meth) %>%
  mutate(days = as_factor(days))

PCA_plot <- ggplot(data = tidy_pca, mapping = aes(x = PC1, y = PC2, label = ID, colour = treatment)) +
  geom_label()


PCA_plot + aes(size = glucose)
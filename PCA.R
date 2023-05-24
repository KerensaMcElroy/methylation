# Load libraries
library(conflicted)
library(tidyverse)
library(readxl)
library(mixOmics)
library(janitor)

# Function to perform PCA analysis and save results
perform_pca <- function(chromosome) {
  # select data for given chromosome
  query <- paste("SELECT methylation, ID, treatment, winMin FROM meth_data WHERE chromosome =", chromosome)
  
  pca_data <- dbGetQuery(con, query) %>%
    as_tibble()
  
  pca_wide <- pca_data %>% pivot_wider(values_from = "methylation", names_from = "winMin")
  pca_wide <- pca_wide %>% 
    dplyr::select(where(~mean(is.na(.)) < .2))
  
  var_df <- pca_wide %>% 
    dplyr::select(-ID,-treatment) %>%
    summarise(across(.cols = everything(), ~var(.x, na.rm = TRUE))) %>% 
    pivot_longer(cols = everything()) %>% 
    arrange(desc(value)) %>%
    mutate(order = 1:n())
  
  var_5000 <- var_df %>%
    slice_max(n = 5000, order_by = value)
  
  hist_plot <- var_df %>% 
    mutate(included = case_when(order <= 5000 ~ TRUE, .default = FALSE)) %>%
    ggplot(aes(x = value, fill = included)) +
    geom_histogram(bins = 100)  +
    scale_fill_brewer(palette = "Accent")
  
  hist_plot
  
  for_pca <- pca_wide %>% 
    dplyr::select(var_5000$name)
  
  result.pca.multi <- for_pca %>%
    pca(scale = TRUE, ncomp = 6, center = TRUE)
  
  # Save PCA values in tidy format
  tidy_pca <- as_tibble(result.pca.multi[["x"]]) %>%
    add_column(ID = pca_wide$ID, treatment = pca_wide$treatment) 
 
  
  # Save the tidy_pca table as a CSV file
  tidy_pca_file <- paste0("results/tidy_pca_", chromosome, ".csv")
  write_csv(tidy_pca, tidy_pca_file)
  
  # Save explained variance
  explained_variance <- as_tibble(result.pca.multi[["explained_variance"]]) %>%
    rownames_to_column(var = "Component")  
  # Save the explained_variance table as a CSV file
  explained_variance_file <- paste0("results/explained_variance_", chromosome, ".csv")
  write_csv(explained_variance, explained_variance_file)
}

# connect to methylation database
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_means/meth_db.sqlite")

# Create a directory to store the results
dir.create("results", showWarnings = FALSE)

# Perform analysis for each chromosome
for (chromosome in 1:29) {
  perform_pca(chromosome)
}

# Disconnect from the database
dbDisconnect(con)


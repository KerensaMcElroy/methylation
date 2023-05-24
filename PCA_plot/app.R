library(shiny)
library(DBI)
library(tidyverse)
library(readxl)
library(mixOmics)
library(janitor)

# UI
ui <- fluidPage(
  titlePanel("Chromosome Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      numericInput("chromosome", "Chromosome:", value = 29, min = 1, max = 29),
      actionButton("runAnalysisBtn", "Run Analysis")
    ),
    mainPanel(
      plotOutput("histPlot"),
      sidebarPanel(
        selectInput("xAxis", "X-Axis:", choices = c("PC1", "PC2", "PC3", "PC4"), selected = "PC1"),
        selectInput("yAxis", "Y-Axis:", choices = c("PC1", "PC2", "PC3", "PC4"), selected = "PC2"),
        selectInput("colorFeature", "Color Feature:", choices = c(), selected = NULL)
      ),
      plotOutput("pcaPlot")
    )
  )
)

# Server
server <- function(input, output) {
  # Reactive value for storing tidy_pca
  tidy_pca_data <- reactiveVal()
  
  # Perform the analysis for a given chromosome
  analyzeData <- function(chromosome) {
    # Path to data
    data_path <- "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/biochem/CC7+8  biochem for Yutao 05_02_21.xlsx"
    
    # Load features data
    features <- data_path %>%
      read_excel(col_names = FALSE) %>% 
      dplyr::select(-1:-2) %>% 
      slice(1:3) %>%
      t() %>%
      as_tibble() %>%
      rename('ID' = 1, 'cohort' = 2, 'room' = 3) %>%
      mutate_if(is.numeric, as.character)
    
    # Load biochem data
    col_names <- read_excel(path = data_path, sheet = 2, range = "A1:Y1", col_names = TRUE)
    biochem <- data_path %>%
      excel_sheets() %>%
      set_names() %>% 
      keep(names(.) != "cheat sheet") %>%
      map_df(~ read_excel(path = data_path, sheet = .x, range =  anchored("A4", c(13, 25)), na='.', col_names = names(col_names)), .id = "biochem") %>%
      rename('days' = 'ID') %>%
      mutate(treatment = case_when(
        days > 17 ~ "pens",
        days > 12 ~ "recovery",
        days > 5  ~ "hot",
        days <= 5 ~ "pre"
      )) %>%
      pivot_longer(cols = 3:26, names_to = 'ID') %>%
      full_join(features, by = 'ID') %>%
      pivot_wider(names_from = 'biochem', values_from = 'value')
    
    # connect to methylation database
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/DIL041/window_read_means/meth_db.sqlite")
    # select data for given chromosome
    query <- paste0("SELECT methylation, ID, treatment, winMin FROM your_table WHERE chromosome =", chromosome)
    
    pca_data <- dbGetQuery(con, query) %>%
      as_tibble()
    
    dbDisconnect(con)
    
    pca_wide <- pca_data %>% pivot_wider(values_from = "methylation", names_from = "winMin")
    pca_wide <- pca_wide %>% 
      dplyr::select(where(~ mean(is.na(.)) < 0.2))
    
    var_df <- pca_wide %>% 
      dplyr::select(-ID, -treatment) %>%
      summarise(across(.cols = everything(), ~ var(.x, na.rm = TRUE))) %>% 
      pivot_longer(cols = everything()) %>% 
      arrange(desc(value)) %>%
      mutate(order = 1:n())
    
    var_1000 <- var_df %>%
      slice_max(n = 1000, order_by = value)
    
    for_pca <- pca_wide %>% 
      dplyr::select(var_1000$name)
    tidy_pca <- as_tibble(result.pca.multi[["x"]]) %>%
      add_column(ID = pca_wide$ID, treatment = pca_wide$treatment) %>%
      inner_join(biochem) %>%
      mutate(days = as_factor(days))
    
    tidy_pca_data(tidy_pca)  # Update the reactive value
  }
  
  # Render histogram plot
  output$histPlot <- renderPlot({
    hist_plot <- var_df %>% 
      mutate(included = case_when(order <= 1000 ~ TRUE, .default = FALSE)) %>%
      ggplot(aes(x = value, fill = included)) +
      geom_histogram(bins = 100) +
      scale_fill_brewer(palette = "Accent")
    
    hist_plot
  })
  
  # Update axis and color feature choices when tidy_pca is available
  observeEvent(tidy_pca_data(), {
    output$xAxis <- renderUI({
      selectInput("xAxis", "X-Axis:", choices = c("PC1", "PC2", "PC3", "PC4"), selected = "PC1")
    })
    
    output$yAxis <- renderUI({
      selectInput("yAxis", "Y-Axis:", choices = c("PC1", "PC2", "PC3", "PC4"), selected = "PC2")
    })
    
    output$colorFeature <- renderUI({
      colorFeatures <- setdiff(colnames(tidy_pca_data()), c("PC1", "PC2", "PC3", "PC4"))
      selectInput("colorFeature", "Color Feature:", choices = colorFeatures, selected = NULL)
    })
  })
  
  # Perform analysis when "Run Analysis" button is clicked
  observeEvent(input$runAnalysisBtn, {
    analyzeData(input$chromosome)
  })
  
  # Render PCA plot
  output$pcaPlot <- renderPlot({
    tidy_pca <- tidy_pca_data()
    if (!is.null(tidy_pca)) {
      PCA_plot <- ggplot(data = tidy_pca, mapping = aes_string(x = input$xAxis, y = input$yAxis, label = "ID", colour = input$colorFeature)) +
        geom_label()
      
      PCA_plot
    }
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

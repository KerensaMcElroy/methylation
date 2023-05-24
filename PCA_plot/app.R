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
      numericInput("chromosome", "Chromosome:", value = 29, min = 1, max = 99),
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
    
   

^    tidy_pca <- as_tibble(result.pca.multi[["x"]]) %>%
      add_column(ID = pca_wide$ID, treatment = pca_wide$treatment) %>%
      inner_join(biochem) %>%
      mutate(days = as_factor(days))
    
    tidy_pca_data(tidy_pca)  # Update the reactive value
  }
  

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

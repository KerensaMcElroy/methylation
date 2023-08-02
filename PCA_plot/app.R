library(shiny)
library(tidyverse)
library(readxl)
library(janitor)


# UI
ui <- fluidPage(
  titlePanel("Chromosome Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      numericInput("chromosome", "Chromosome:", value = 1, min = 1, max = 29),
      actionButton("runAnalysisBtn", "Run Analysis"),
      br(),
      br(),
      selectInput("xAxis", "X-Axis:", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), selected = "PC1"),
      selectInput("yAxis", "Y-Axis:", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), selected = "PC2"),
      selectInput("colorFeature", "Color Feature:", choices = c("treatment", "days", "cohort",
                                                                "room", "sodium", "potassium", "chloride", "bicarb", "anion_gap", "Mg",
                                                                "glucose", "BOHB", "cholesterol", "urea", "creatinine", "U2C_ratio", "Ca",
                                                                "Phos", "Ca2P_ratio", "T_protein", "albumin", "globulins", "A2G_ratio",
                                                                "bilirubin", "ALP", "AST", "CK", "GTT", "GLDH", "glutamine", "glutamate"), selected = "treatment"
    ),
      selectInput("windowRange", "Window Range:", 
                choices = c("1-5000", "5001-10000", "10001-15000", "15001-20000"),
                selected = "1-5000")
    ),
  
    
    mainPanel(
      plotOutput("pcaPlot", height = "80vh")
    )
  )
)

# Server
server <- function(input, output) {
  # Reactive values
  tidy_pca_data <- reactiveVal()
  pca_plot <- reactiveVal()
  color_features <- reactiveVal()
  explained_variance <- reactiveVal()
  
  mapChoiceToWindowRange <- function(choice) {
    switch(choice,
           "1-5000"      = ".csv",
           "5001-10000"  = "_5001_10000.csv",
           "10001-15000" = "_10001_15000.csv",
           "15001-20000" = "_15001-20000.csv"
    )
  }
  
  # Perform the analysis for a given chromosome and window range
  analyzeData <- function(chromosome, window_choice) {
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
      pivot_wider(names_from = 'biochem', values_from = 'value') %>%
      dplyr::filter(days %in% c(5, 12, 38))  # to match meth collection events
    
    window_range <- mapChoiceToWindowRange(window_choice)

    pca_data <- read_csv(paste0("results/tidy_pca_", chromosome, window_range))
    pca_data <- pca_data %>% mutate(ID = as.character(ID))
    
    explained_variance_data <- read_csv(paste0("results/explained_variance_", chromosome, window_range))
    
    explained_variance(explained_variance_data)
    
    # Print structure of pca_data
    print("Structure of pca_data:")
    str(pca_data)
    
    # Print structure of biochem
    print("Structure of biochem:")
    str(biochem)
    
    tidy_pca <- pca_data %>%
      inner_join(biochem, by = c("ID", "treatment")) %>%
      mutate(days = as_factor(days))
    
    print("Structure of tidy:")
    str(tidy_pca)
    
    tidy_pca_data(tidy_pca)
    pca_plot(NULL)
    
  }
  
  
  # Update axis and color feature choices when tidy_pca is available
  observeEvent(tidy_pca_data(), {
    output$xAxis <- renderUI({
      selectInput("xAxis", "X-Axis:", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), selected = "PC1")
    })
    
    output$yAxis <- renderUI({
      selectInput("yAxis", "Y-Axis:", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), selected = "PC2")
    })
    
    output$colorFeature <- renderUI({
      choices <- c("treatment", "days", "cohort",
                   "room", "sodium", "potassium", "chloride", "bicarb", "anion_gap", "Mg",
                   "glucose", "BOHB", "cholesterol", "urea", "creatinine", "U2C_ratio", "Ca",
                   "Phos", "Ca2P_ratio", "T_protein", "albumin", "globulins", "A2G_ratio",
                   "bilirubin", "ALP", "AST", "CK", "GTT", "GLDH", "glutamine", "glutamate")
      selectInput("colorFeature", "Color Feature:", choices = choices, selected = "treatment")
    })
    
  })
  # Perform analysis when "Run Analysis" button is clicked
  observeEvent(input$runAnalysisBtn, {
    analyzeData(input$chromosome, input$windowRange)
  })
  observeEvent(input$chromosome, {
    analyzeData(input$chromosome, input$windowRange)
  })
  observeEvent(input$windowRange, {
    analyzeData(input$chromosome, input$windowRange)
  })
  # Render PCA plot
  output$pcaPlot <- renderPlot({
    tidy_pca <- tidy_pca_data()
    if (!is.null(tidy_pca)) {
      x_axis <- input$xAxis
      y_axis <- input$yAxis
      
      pca_plot(NULL)  # Clear the previous plot
      
      plot_title <- paste0("PCA Analysis: ", x_axis, " vs ", y_axis)
      
      x_variance <- round(explained_variance()$value[explained_variance()$Component == as.numeric(str_extract(x_axis, "\\d+"))] * 100, 2)
      y_variance <- round(explained_variance()$value[explained_variance()$Component == as.numeric(str_extract(y_axis, "\\d+"))] * 100, 2)
      x_label <- paste0(x_axis, " (Explained Variance: ", x_variance, "%)")
      y_label <- paste0(y_axis, " (Explained Variance: ", y_variance, "%)")
      
      color_feature <- input$colorFeature
      pca_plot(
        ggplot(data = tidy_pca, mapping = aes_string(x = x_axis, y = y_axis, label = "ID", fill = color_feature)) +
          geom_label(color="black") +
          labs(title = plot_title, x = x_label, y = y_label)
      )
      
      pca_plot()
    }
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

# app.R
library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggridges)
library(reshape2)
library(dplyr)
library(ranger)
library(pROC)
library(rstatix)
library(purrr)
library(rsconnect)
library(patchwork)

# Source utility functions
source("utils.R")

# Default data (loaded only if no upload)
shiny_data <- readRDS("shiny_data.Rds")
shiny_meta <- readRDS("shiny_meta.Rds")

# UI
ui <- fluidPage(
  title = "MachOmics",
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
  tags$div(class = "navbar navbar-fixed-top",
           tags$a(href = "https://pdobrano25.github.io/ml_website/", 
                  tags$img(src = "machomics_flat_white.png", alt = "MachOmics", style = "height:30px;")),
           tags$ul(class = "nav navbar-nav navbar-right",
                   tags$li(tags$a(href = "https://pdobrano25.github.io/ml_website/", "Home")),
                   tags$li(tags$a(href = "https://pdobrano25.github.io/ml_website/machomics.html", "MachOmics", class = "active")),
                   tags$li(tags$a(href = "https://pdobrano25.github.io/ml_website/about.html", "About"))
           )
  ),
  fluidRow(
    column(12, h1("MachOmics"))
  ),
  # Data Upload Section
  fluidRow(
    column(12, h4("Data Upload Section", style = "font-weight: bold; margin-bottom: 10px;"),
           p("Upload your feature data and metadata in CSV, TSV, or RDS format, or use the default datasets. Configure the task type and labels, and apply transformations to prepare your data for analysis. You can also download the default data as CSV files below.")),
    box(
      title = "Data Upload", width = 6, status = "primary",
      fileInput("data_upload", "Upload Feature Data (CSV/TSV/RDS)", accept = c(".csv", ".tsv", ".rds", ".Rds", ".RDS"), placeholder = "Drag and drop or click to upload"),
      fileInput("meta_upload", "Upload Metadata (CSV/TSV/RDS)", accept = c(".csv", ".tsv", ".rds", ".Rds", ".RDS"), placeholder = "Drag and drop or click to upload"),
      downloadButton("download_default_data", "Download Default Feature Data (CSV)"),
      downloadButton("download_default_meta", "Download Default Metadata (CSV)"),
      verbatimTextOutput("upload_status")
    ),
    box(
      title = "Data Configuration", width = 6, status = "primary",
      selectInput("data_type", "Task Type", choices = c("classification", "regression")),
      textInput("case_label", "Case Label", value = "CD"),
      textInput("control_label", "Control Label", value = "no"),
      actionButton("update_meta", "Update Metadata"),
      checkboxInput("apply_scale", "Scale to 100%", value = FALSE),
      checkboxInput("apply_transform", "Apply Log2 Transform", value = TRUE),
      actionButton("update_transform", "Apply Transformation")
    )
  ),
  # Data Exploration Section
  fluidRow(
    column(12, h4("Data Exploration Section", style = "font-weight: bold; margin-bottom: 10px;"),
           p("Visualize the distribution of your original and transformed feature data to understand its characteristics. Check the target variable distribution to assess class imbalances or regression ranges. Download the plots for further analysis.")),
    box(
      title = "Original Distribution", width = 6, status = "info",
      plotOutput("dist_plot_original"),
      downloadButton("download_dist_original", "Download Original Plot")
    ),
    box(
      title = "Transformed Distribution", width = 6, status = "info",
      plotOutput("dist_plot_transformed"),
      downloadButton("download_dist_transformed", "Download Transformed Plot")
    )
  ),
  fluidRow(
    box(
      title = "Target Distributions", width = 12, status = "info",
      plotOutput("target_plot"),
      downloadButton("download_target", "Download Plot")
    )
  ),
  # Model Building Section
  fluidRow(
    column(12, h4("Model Building Section", style = "font-weight: bold; margin-bottom: 10px;"),
           p("Configure and train a machine learning model for classification or regression using your processed data. Monitor the training progress and status to ensure successful model building.")),
    box(
      title = "Model Configuration", width = 6, status = "warning",
      selectInput("model_task", "Task Type", choices = c("classification", "regression")),
      actionButton("run_model", "Run Model")
    ),
    box(
      title = "Model Status", width = 6, status = "warning",
      verbatimTextOutput("model_status")
    )
  ),
  # Results Section
  fluidRow(
    column(12, h4("Results Section", style = "font-weight: bold; margin-bottom: 10px;"),
           p("Review the model’s performance metrics and feature importances through interactive plots. Download these visualizations to share or include in reports.")),
    box(
      title = "Performance Metrics", width = 12, status = "success",
      plotOutput("perf_plot", width = "1200px", height = "400px"),
      downloadButton("download_perf", "Download Plot")
    )
  ),
  fluidRow(
    box(
      title = "Feature Importances", width = 12, status = "success",
      plotOutput("feat_plot", width = "1200px", height = "400px"),
      downloadButton("download_feat", "Download Plot")
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive values to store processed data
  rv <- reactiveValues(
    meta = shiny_meta,
    data = shiny_data,
    data_original = shiny_data,  # Store original untransformed data
    model_output = NULL,
    upload_status = "Using default data",
    transformed = FALSE  # Track if transformation has been applied
  )
  
  # Handle Data Upload
  observeEvent(input$data_upload, {
    req(input$data_upload)
    tryCatch({
      file_ext <- tools::file_ext(input$data_upload$datapath)
      if (file_ext == "csv") {
        rv$data <- read.csv(input$data_upload$datapath, row.names = 1, check.names = FALSE)
      } else if (file_ext %in% c("tsv", "txt")) {
        rv$data <- read.delim(input$data_upload$datapath, row.names = 1, check.names = FALSE)
      } else if (file_ext %in% c("rds", "Rds", "RDS")) {
        rv$data <- readRDS(input$data_upload$datapath)
      } else {
        stop("Unsupported file format. Please upload CSV, TSV, or RDS.")
      }
      # Validate data structure
      if (!is.matrix(rv$data) && !is.data.frame(rv$data)) {
        stop("Uploaded data must be a matrix or dataframe with samples as rows.")
      }
      if (!all(sapply(rv$data, is.numeric))) {
        stop("Feature data must contain only numeric values.")
      }
      # Store original data and reset transformation
      rv$data_original <- rv$data
      rv$transformed <- FALSE
      rv$upload_status <- paste("Feature data uploaded successfully:", input$data_upload$name)
    }, error = function(e) {
      rv$upload_status <- paste("Error uploading feature data:", e$message)
    })
  })
  
  # Handle Meta Upload
  observeEvent(input$meta_upload, {
    req(input$meta_upload)
    tryCatch({
      file_ext <- tools::file_ext(input$meta_upload$datapath)
      if (file_ext == "csv") {
        rv$meta <- read.csv(input$meta_upload$datapath, check.names = FALSE)
      } else if (file_ext == "tsv") {
        rv$meta <- read.delim(input$meta_upload$datapath, check.names = FALSE)
      } else if (file_ext %in% c("rds", "Rds", "RDS")) {
        rv$meta <- readRDS(input$meta_upload$datapath)
      } else {
        stop("Unsupported file format. Please upload CSV, TSV, or RDS.")
      }
      # Validate meta structure
      if (!all(c("sample", "value") %in% colnames(rv$meta))) {
        stop("Metadata must have 'sample' and 'value' columns.")
      }
      rv$upload_status <- paste(rv$upload_status, "\nMetadata uploaded successfully:", input$meta_upload$name)
    }, error = function(e) {
      rv$upload_status <- paste("Error uploading metadata:", e$message)
    })
  })
  
  # Download Default Feature Data as CSV
  output$download_default_data <- downloadHandler(
    filename = function() { "default_feature_data.csv" },
    content = function(file) {
      write.csv(shiny_data, file, row.names = TRUE)
    }
  )
  
  # Download Default Metadata as CSV
  output$download_default_meta <- downloadHandler(
    filename = function() { "default_metadata.csv" },
    content = function(file) {
      write.csv(shiny_meta, file, row.names = FALSE)
    }
  )
  
  # Display Upload Status
  output$upload_status <- renderText({
    rv$upload_status
  })
  
  # Update Metadata
  observeEvent(input$update_meta, {
    tryCatch({
      rv$meta <- change_names(
        meta = rv$meta,
        type = input$data_type,
        case = input$case_label,
        control = input$control_label
      )
      rv$upload_status <- paste(rv$upload_status, "\nMetadata updated successfully")
    }, error = function(e) {
      rv$upload_status <- paste("Error updating metadata:", e$message)
    })
  })
  
  # Update Transformation
  observeEvent(input$update_transform, {
    tryCatch({
      rv$data <- log_transform(
        data = rv$data_original,
        scale = input$apply_scale,
        transform = input$apply_transform
      )
      rv$transformed <- input$apply_transform || input$apply_scale
      rv$upload_status <- paste(rv$upload_status, "\nData transformation applied successfully")
    }, error = function(e) {
      rv$upload_status <- paste("Error applying transformation:", e$message)
    })
  })
  
  # Original Distribution Plot
  output$dist_plot_original <- renderPlot({
    req(rv$data_original)
    distribution_check(data = rv$data_original, transform = FALSE)
  })
  
  # Download Original Distribution Plot
  output$download_dist_original <- downloadHandler(
    filename = function() { "distribution_plot_original.png" },
    content = function(file) {
      ggsave(file, plot = distribution_check(data = rv$data_original, transform = FALSE), device = "png")
    }
  )
  
  # Transformed Distribution Plot
  output$dist_plot_transformed <- renderPlot({
    req(rv$transformed)
    distribution_check(data = rv$data_original, transform = input$apply_transform)
  })
  
  # Download Transformed Distribution Plot
  output$download_dist_transformed <- downloadHandler(
    filename = function() { "distribution_plot_transformed.png" },
    content = function(file) {
      req(rv$transformed)
      ggsave(file, plot = distribution_check(data = rv$data_original, transform = input$apply_transform), device = "png")
    }
  )
  
  # Target Plot
  output$target_plot <- renderPlot({
    check_imbalances(meta = rv$meta, type = input$data_type)
  })
  
  # Download Target Plot
  output$download_target <- downloadHandler(
    filename = function() { "target_plot.png" },
    content = function(file) {
      ggsave(file, plot = check_imbalances(meta = rv$meta, type = input$data_type), device = "png")
    }
  )
  
  # Run Model with Progress Bar
  observeEvent(input$run_model, {
    tryCatch({
      withProgress(message = "Training model...", value = 0, {
        rv$model_output <- build_model(
          data = rv$data,
          meta = rv$meta,
          task = input$model_task
        )
        incProgress(1)
      })
      output$model_status <- renderText({"Model training completed successfully!"})
    }, error = function(e) {
      output$model_status <- renderText({
        paste("Error training model:", e$message)
      })
    })
  })
  
  # Performance Plot
  output$perf_plot <- renderPlot({
    req(rv$model_output)
    visualize_performance(build_output = rv$model_output, task = input$model_task)
  })
  
  # Download Performance Plot
  output$download_perf <- downloadHandler(
    filename = function() { "performance_plot.png" },
    content = function(file) {
      req(rv$model_output)
      ggsave(file, plot = visualize_performance(build_output = rv$model_output, task = input$model_task), 
             device = "png", width = 1200/96, height = 400/96, units = "in")
    }
  )
  
  # Feature Importance Plot
  output$feat_plot <- renderPlot({
    req(rv$model_output)
    withProgress(message = "Generating feature importances...", value = 0, {
      plot <- feature_importances(
        build_output = rv$model_output,
        data = rv$data,
        meta = rv$meta
      )
      incProgress(1)
      plot
    })
  )
  
  # Download Feature Importance Plot
  output$download_feat <- downloadHandler(
    filename = function() { "feature_importance_plot.png" },
    content = function(file) {
      req(rv$model_output)
      ggsave(file, plot = feature_importances(build_output = rv$model_output, data = rv$data, meta = rv$meta), 
             device = "png", width = 1200/96, height = 400/96, units = "in")
    }
  )
  }
  
  # Run App
  shinyApp(ui, server)
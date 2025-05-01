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
library(bslib)  # For Cosmo theme
library(htmltools)  # For includeHTML

# deploy
# setwd("~/Documents/PhD/ml_website")

# rsconnect::setAccountInfo(name = "pdobrano", token = "DEDBB05D0EA0E7BC39EEB71CC0FAE1FB", secret = "/6KZAFGXuzLFXcRpAy/KcZ4b+MfMWgI1wec5OhCo")

# rsconnect::deployApp('MachOmics')


# setwd("MachOmics") # run this before testing
# shiny::runApp("app.R") # use this to test before deploying

source("utils.R")

# Default data (loaded only if no upload)
shiny_data <- readRDS("shiny_data.Rds")
shiny_meta <- readRDS("shiny_meta.Rds")

# :: navbar ----------------------------------------------------------------

# add root
addResourcePath("root", ".")

navbar_html <-({
  html_content <- readLines("_includes/navbar.html")
  html_content <- gsub('src\\s*=\\s*[\'"]machomics_flat_white\\.png[\'"]', 'src="root/machomics_flat_white.png"', html_content, ignore.case = TRUE)
  html_content <- gsub('src="www/machomics_flat_white\\.png"', 'src="root/machomics_flat_white.png"', html_content, ignore.case = TRUE)
  if (!file.exists("machomics_flat_white.png")) {
    cat("Root image missing. Using www/ or fallback URL\n")
    html_content <- gsub('src="root/machomics_flat_white\\.png"', 'src="www/machomics_flat_white.png"', html_content, ignore.case = TRUE)
    if (!file.exists("www/machomics_flat_white.png")) {
      cat("www/ image missing. Using GitHub URL\n")
      html_content <- gsub('src="www/machomics_flat_white\\.png"', 'src="https://pdobrano25.github.io/ml_website/images/machomics_flat_white.png"', html_content, ignore.case = TRUE)
    } else {
      cat("Using www image: www/machomics_flat_white.png\n")
    }
  } else {
    cat("Using root image: root/machomics_flat_white.png\n")
  }
  html_content <- gsub('href=\\"(?!(http|https|ftp|#))([^"]*)"', 'href="https://pdobrano25.github.io/ml_website/\\2"', html_content, perl = TRUE)
  paste(html_content, collapse="\n")
})

# :: UI -------------------------------------------------------------------

ui <- fluidPage(
  theme = bs_theme(bootswatch = "cosmo"),
  title = "MachOmics",
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$style(HTML("
      body { 
        padding-top: 70px;
        padding-left: 220px;
      }
      .navbar { z-index: 1000; }
      .navbar-brand img { 
        display: inline-block !important; 
        visibility: visible !important; 
        height: 30px !important; 
        max-width: 150px;
        object-fit: contain;
      }
      .download-button { opacity: 1 !important; pointer-events: auto !important; }
      .navbar-default {
        background-color: #1A1A1A !important;
        border-color: #222222 !important;
      }
      .navbar-default .navbar-brand,
      .navbar-default .navbar-nav > li > a {
        color: #FFFFFF !important;
      }
      .navbar-default .navbar-nav > li > a:hover,
      .navbar-default .navbar-nav > li.active > a {
        color: #1A1A1A !important;
        background-color: #FFFFFF !important;
      }
      .navbar-default .dropdown-menu {
        background-color: #1A1A1A !important;
      }
      .navbar-default .dropdown-menu > li > a {
        color: #FFFFFF !important;
      }
      .navbar-default .dropdown-menu > li > a:hover,
      .navbar-default .dropdown-menu > li.active > a {
        color: #1A1A1A !important;
        background-color: #FFFFFF !important;
      }
      .navbar-default .navbar-toggle {
        border-color: transparent !important;
        background-color: transparent !important;
      }
      .navbar-default .navbar-toggle .icon-bar {
        background-color: #FFFFFF !important;
      }
      .navbar-default .navbar-toggle:hover,
      .navbar-default .navbar-toggle:focus {
        background-color: #555555 !important;
      }
      #toc {
        position: fixed;
        top: 70px;
        left: 0;
        width: 200px;
        background-color: #FFFFFF;
        border-right: 1px solid #DDD;
        padding: 10px;
        max-height: calc(100vh - 70px);
        overflow-y: auto;
        z-index: 999;
      }
      #toc ul {
        list-style: none;
        padding-left: 10px;
      }
      #toc li {
        margin-bottom: 5px;
      }
      #toc a {
        color: #000000;
        background-color: #FFFFFF;
        text-decoration: none;
        font-size: 14px;
        font-weight: bold;
        display: block;
        padding: 5px 10px;
        border-radius: 4px;
      }
      #toc a:hover {
        color: #000000;
        background-color: #F5F5F5;
      }
      #toc a.active {
        color: #FFFFFF;
        background-color: #337AB7;
      }
      @media (max-width: 768px) {
        body { padding-left: 0; }
        #toc { display: none; }
      }
    ")),
    tags$script(HTML("
      document.addEventListener('DOMContentLoaded', function() {
        const tocLinks = document.querySelectorAll('#toc a');
        const sections = document.querySelectorAll('h4[id]');
        tocLinks.forEach(link => {
          link.addEventListener('click', function() {
            tocLinks.forEach(l => l.classList.remove('active'));
            this.classList.add('active');
          });
        });
        window.addEventListener('scroll', function() {
          let current = '';
          sections.forEach(section => {
            const sectionTop = section.offsetTop;
            if (window.pageYOffset >= sectionTop - 70) {
              current = section.getAttribute('id');
            }
          });
          tocLinks.forEach(link => {
            link.classList.remove('active');
            if (link.getAttribute('href') === '#' + current) {
              link.classList.add('active');
            }
          });
        });
      });
    "))
  ),
  HTML(navbar_html),
  tags$div(id = "toc",
           tags$ul(
             tags$li(tags$a(href = "#data-upload-section", "Data Upload Section")),
             tags$li(tags$a(href = "#data-exploration-section", "Data Exploration Section")),
             tags$li(tags$a(href = "#model-building-section", "Model Building Section")),
             tags$li(tags$a(href = "#results-section", "Results Section"))
           )
  ),
  fluidRow(
    column(12, h1("MachOmics"))
  ),
  # Data Upload Section
  fluidRow(
    column(12, h4("Data Upload Section", style = "font-weight: bold; margin-bottom: 10px;", id = "data-upload-section"),
           p("Upload your feature data and metadata in CSV, TSV, or RDS format, or use the default datasets. Configure the task type and labels, and apply transformations to prepare your data for analysis. You can also download the default data as CSV files below.")),
    box(
      title = "Data Upload", width = 6, status = "primary",
      fileInput("data_upload", "Upload Feature Data (CSV/TSV/RDS)", accept = c(".csv", ".tsv", ".rds", ".Rds", ".RDS"), placeholder = "...or drag and drop"),
      fileInput("meta_upload", "Upload Metadata (CSV/TSV/RDS)", accept = c(".csv", ".tsv", ".rds", ".Rds", ".RDS"), placeholder = "...or drag and drop"),
      downloadButton("download_default_data", "Download Default Feature Data (CSV)", style = "margin-top: 10px;", class = "download-button"),
      downloadButton("download_default_meta", "Download Default Metadata (CSV)", style = "margin-top: 10px;", class = "download-button"),
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
      numericInput("prevalence", "Filter Features: Prevalence Threshold (%)", value = 10, min = 0, max = 100, step = 1),
      actionButton("apply_filter", "Apply Feature Filter"),
      actionButton("update_transform", "Apply Transformation"),
      verbatimTextOutput("filter_status")
    )
  ),
  # Data Exploration Section
  fluidRow(
    column(12, h4("Data Exploration Section", style = "font-weight: bold; margin-bottom: 10px;", id = "data-exploration-section"),
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
      title = "Target Distributions", width = 6, status = "info",
      plotOutput("target_plot"),
      downloadButton("download_target", "Download Plot")
    )
  ),
  # Model Building Section
  fluidRow(
    column(12, h4("Model Building Section", style = "font-weight: bold; margin-bottom: 10px;", id = "model-building-section"),
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
    column(12, h4("Results Section", style = "font-weight: bold; margin-bottom: 10px;", id = "results-section"),
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

# :: Server -------------------------------------------------------------------

server <- function(input, output, session) {
  # Reactive values to store processed data
  rv <- reactiveValues(
    meta = shiny_meta,
    data = shiny_data,
    data_original = shiny_data,  # Store original untransformed data
    model_output = NULL,
    upload_status = "Using default data",
    transformed = FALSE,  # Track if transformation has been applied
    filtered = FALSE,  # Track if filter has been applied
    feature_count_before = ncol(shiny_data),  # Initial feature count
    feature_count_after = ncol(shiny_data)  # Feature count after filtering
  )
  
  # Download sample data
  output$download_default_data <- downloadHandler(
    filename = function() { "default_feature_data.csv" },
    content = function(file) {
      write.csv(shiny_data, file, row.names = TRUE)
    }
  )
  
  # Download sample metadata
  output$download_default_meta <- downloadHandler(
    filename = function() { "default_metadata.csv" },
    content = function(file) {
      write.csv(shiny_meta, file, row.names = FALSE)
    }
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
      if (!is.matrix(rv$data) && !is.data.frame(rv$data)) {
        stop("Uploaded data must be a matrix or dataframe with samples as rows.")
      }
      if (!all(sapply(rv$data, is.numeric))) {
        stop("Feature data must contain only numeric values.")
      }
      # Store original data and reset transformation/filter
      rv$data_original <- rv$data
      rv$transformed <- FALSE
      rv$filtered <- FALSE
      rv$feature_count_before <- ncol(rv$data)
      rv$feature_count_after <- ncol(rv$data)
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
      if (!all(c("sample", "value") %in% colnames(rv$meta))) {
        stop("Metadata must have 'sample' and 'value' columns.")
      }
      rv$upload_status <- paste(rv$upload_status, "\nMetadata uploaded successfully:", input$meta_upload$name)
    }, error = function(e) {
      rv$upload_status <- paste("Error uploading metadata:", e$message)
    })
  })
  
  # Display Upload Status
  output$upload_status <- renderText({
    rv$upload_status
  })
  
  # Apply Feature Filter
  observeEvent(input$apply_filter, {
    tryCatch({
      req(rv$data_original)
      rv$feature_count_before <- ncol(rv$data_original)
      rv$data <- filter_data(data = rv$data_original, prev = input$prevalence)
      rv$feature_count_after <- ncol(rv$data)
      rv$filtered <- TRUE
      rv$upload_status <- paste(
        rv$upload_status,
        sprintf(
          "\nFeature filter applied (prevalence: %d%%). Features before: %d, after: %d",
          input$prevalence, rv$feature_count_before, rv$feature_count_after
        )
      )
    }, error = function(e) {
      rv$upload_status <- paste("Error applying feature filter:", e$message)
    })
  })
  
  # Display Filter Status
  output$filter_status <- renderText({
    if (rv$filtered) {
      sprintf(
        "Features before filtering: %d\nFeatures after filtering: %d (prevalence: %d%%)",
        rv$feature_count_before, rv$feature_count_after, input$prevalence
      )
    } else {
      sprintf("No filter applied. Total features: %d", rv$feature_count_before)
    }
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
      # Apply transformation to filtered data (rv$data) if filtered, else original
      data_to_transform <- if (rv$filtered) rv$data else rv$data_original
      rv$data <- log_transform(
        data = data_to_transform,
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
      ggsave(file, plot = distribution_check(data = rv$data_original, transform = FALSE), 
             width = 6, height = 4, units = "in", dpi = 320, device = "png")
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
      ggsave(file, plot = distribution_check(data = rv$data_original, transform = input$apply_transform), 
             width = 6, height = 4, units = "in", dpi = 320, device = "png")
    }
  )
  
  # Target Plot
  output$target_plot <- renderPlot({
    req(rv$meta)
    check_imbalances(meta = rv$meta, type = input$data_type)
  })
  
  # Download Target Plot
  output$download_target <- downloadHandler(
    filename = function() { "target_plot.png" },
    content = function(file) {
      ggsave(file, plot = check_imbalances(meta = rv$meta, type = input$data_type), 
             width = 6, height = 4, units = "in", dpi = 320, device = "png")
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
      ggsave(file, plot = visualize_performance(build_output = rv$model_output, task = input$model_task), 
             width = 12, height = 4, units = "in", dpi = 320, device = "png")
    }
  )
  
  # Feature Importance Plot with Progress Bar
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
  })
  
  # Download Feature Importance Plot
  output$download_feat <- downloadHandler(
    filename = function() { "feature_importance_plot.png" },
    content = function(file) {
      ggsave(file, plot = feature_importances(build_output = rv$model_output, data = rv$data, meta = rv$meta),
             width = 12, height = 4, units = "in", dpi = 320, device = "png")
    }
  )
}

# Run App
shinyApp(ui, server)
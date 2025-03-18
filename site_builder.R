### Constructing ml_website

load("2025_03_18_ml_website_image.Renv")

#.  save.image("2025_03_18_ml_website_image.Renv")

# :: render navbar --------------------------------------------------------

# NOTE: Don't forget to update _site.yml with new pages in navbar

library(yaml)

# Read the _site.yml file
site_config <- yaml::read_yaml("_site.yml")

# Extract navbar
navbar <- site_config$navbar
# Temporary render to get navbar HTML
unlink("docs", recursive = TRUE)
# use index to render navbar (saved to root, not)
rmarkdown::render(input = "index.Rmd",
                  output_options = list(
                         theme = "cosmo"))
html_content <- readLines("index.html", warn = FALSE)
nav_start <- grep("<div.*navbar", html_content, ignore.case = TRUE)[1]
nav_end <- grep("<p>", html_content[nav_start:length(html_content)])[1] + nav_start -2
nav_remove <- grep("<h1 class=\"title toc-ignore\">Machine Learning Projects</h1>", html_content)
nav_index <- c(nav_start:(nav_remove-1), (nav_remove+1):nav_end)
navbar_html <- c(
  html_content[nav_index],
  '<style>',
  '  .navbar-fixed-top { position: fixed; top: 0; width: 100%; z-index: 1000; }',  # Keep fixed behavior
  '  body { padding-top: 70px; }',  # Offset content below navbar (adjust 70px as needed)
  '</style>'
)
# Write navbar.html
dir.create("_includes")
writeLines(navbar_html, "_includes/navbar.html")
navbar.path = normalizePath("_includes/navbar.html", mustWork = TRUE)


# :: render html ---------------------------------------------------------


# List files
files <- c("index.Rmd", "about.Rmd", "mlp_validation/mlp_validation.Rmd",
           "gevers_validation/gevers_validation.Rmd", "ml_models/ml_figures.Rmd")

# Clean output
unlink("docs", recursive = TRUE)

# Render each file
for (f in files) {
  output_dir <- file.path("docs", dirname(f))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Define output options
  output_options <- list(
    theme = "cosmo",
    toc = TRUE,
    toc_float = TRUE
  )
  
  # Conditionally add navbar (otherwise navbar gets doubled)
  if (!(basename(f) %in% c("index.Rmd", "about.Rmd"))) {
    output_options$includes <- list(before_body = navbar.path)
  }
  
  # Render the file
  rmarkdown::render(f, 
                    output_dir = "docs",
                    output_format = "html_document", 
                    output_options = output_options)
}



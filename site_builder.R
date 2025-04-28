### Constructing ml_website

# load("2025_03_18_ml_website_image.Renv")

#.  save.image("2025_03_18_ml_website_image.Renv")

# :: render navbar --------------------------------------------------------

# NOTE: Don't forget to update _site.yml with new pages in navbar

t1 <- Sys.time()
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
nav_remove <- grep("<h1 class=\"title toc-ignore\">MachOmics</h1>", html_content)
nav_index <- c(nav_start:(nav_remove-1), (nav_remove+1):nav_end)

# Modify navbar to replace title with logo
navbar_html <- html_content[nav_index]
title_line <- grep("<a class=\"navbar-brand\" href=\"index.html\">MachOmics</a>", navbar_html)
navbar_html[title_line] <- '<a class="navbar-brand" href="index.html"><img src="machomics_flat_white.png" alt="MachOmics" style="height:30px;"></a>'


navbar_html <- c(
  html_content[nav_index],
  '<style>',
  '  .navbar-fixed-top { position: fixed; top: 0; width: 100%; z-index: 1000; }',  # Keep fixed behavior
  '  body { padding-top: 70px; }',  # Offset content below navbar (adjust 70px as needed)
  '  .navbar-brand img { vertical-align: middle; }',
  '</style>'
)
# Write navbar.html
dir.create("_includes")
writeLines(navbar_html, "_includes/navbar.html")
navbar.path = normalizePath("_includes/navbar.html", mustWork = TRUE)


# :: render html ---------------------------------------------------------


# List files
files <- c("index.Rmd", 
           "MachOmics/machomics.Rmd",
           "about.Rmd", 
           "mlp_validation/mlp_validation.Rmd",
           "gevers_validation/gevers_validation.Rmd", 
           "cmd_validation/cmd_validation.Rmd",
           "sinai_validation/sinai_validation.Rmd",
           "metaaml_validation/metaaml_validation.Rmd",   
           "hmp_validation/hmp_validation.Rmd",
           "muller_validation/muller_validation.Rmd",
           "litichevskiy_validation/litichevskiy_validation.Rmd",
           "overfitcheck.Rmd",
           "ml_models/ml_figures.Rmd")

# Clean output
unlink("docs", recursive = TRUE)

# Render each file
for (f in files) {
  output_dir <- file.path("docs", dirname(f))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Copy images to docs/
  dir.create("docs", recursive = TRUE, showWarnings = FALSE)
  for (img in c("logo-black.png", "linkedin_icon.png", "machomics_flat_white.png")) {
    img_source <- img
    img_dest <- file.path("docs", img)
    if (file.exists(img_source)) {
      file.copy(img_source, img_dest, overwrite = TRUE)
    }
  }
  # Copy logo to MachOmics/ for Shiny app
  dir.create("MachOmics", recursive = TRUE, showWarnings = FALSE)
  if (file.exists("machomics_flat_white.png")) {
    file.copy("machomics_flat_white.png", "MachOmics/machomics_flat_white.png", overwrite = TRUE)
  }
  
  # Define output options
  output_options <- list(
    theme = "cosmo",
    toc = TRUE,
    toc_float = TRUE
  )
  
  # Conditionally add navbar (otherwise navbar gets doubled)
  if (!(basename(f) %in% c("index.Rmd", "overfitcheck.Rmd", "about.Rmd"))) {
    output_options$includes <- list(before_body = navbar.path)
  }
  
  # Render the file
  rmarkdown::render(f, 
                    output_dir = "docs",
                    output_format = "html_document", 
                    output_options = output_options)
}

t2 <- Sys.time()
t2-t1 
# Apr 4 2025 = 12 min (added hmp_validation)
# Apr 10 2025 = 17 min (added muller_validation)
# Apr 12 2025 = 20 min (added litichevskiy_validation)
# Apr 27 2025 = 35 min (added shiny app and logo


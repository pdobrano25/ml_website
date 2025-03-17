### Constructin ml_website


# :: render navbar --------------------------------------------------------

# NOTE: Don't forget to update _site.yml with new pages in navbar

library(yaml)

# Read the _site.yml file
site_config <- yaml::read_yaml("_site.yml")

# Extract navbar
navbar <- site_config$navbar

# Function to generate HTML for a menu item
generate_menu_item <- function(item) {
  if (!is.null(item$menu)) {
    submenu <- paste0(
      '<li class="dropdown">',
      '<a href="#" class="dropbtn">', item$text, ' ▼</a>',
      '<ul class="dropdown-content">',
      paste(sapply(item$menu, function(sub) {
        sprintf('<li><a href="%s">%s</a></li>', sub$href, sub$text)
      }), collapse = ""),
      '</ul>',
      '</li>'
    )
    return(submenu)
  } else {
    return(sprintf('<li><a href="%s">%s</a></li>', item$href, item$text))
  }
}

# Generate HTML navbar with title as a button
navbar_html <- paste0(
  '<nav>',
  '<div class="navbar-header">',
  '<a href="index.html" class="navbar-brand">', navbar$title, '</a>',  # Title as button
  '</div>',
  '<ul>',
  paste(sapply(navbar$left, generate_menu_item), collapse = ""),
  '</ul>',
  '</nav>',
  '<style>
  nav { background: #333; overflow: hidden; padding: 10px; }
  nav ul { list-style-type: none; margin: 0; padding: 0; float: left; }
  nav ul li { float: left; position: relative; }
  nav ul li a { display: block; color: white; padding: 14px 16px; text-decoration: none; }
  nav ul li a:hover { background: #111; }
  
  .navbar-header { float: left; }
  .navbar-brand { display: block; color: white; padding: 14px 16px; text-decoration: none; font-size: 18px; }
  .navbar-brand:hover { background: #111; }
  
  .dropdown-content { 
    display: none; 
    position: absolute; 
    top: 100%; 
    left: 0; 
    background: #333; 
    min-width: 160px; 
    z-index: 1000;
  }
  
  .dropdown-content li { float: none; }
  .dropdown-content li a { padding: 12px 16px; display: block; }
  .dropdown-content li a:hover { background: #555; }
  
  .dropdown:hover .dropdown-content { display: block; }
  </style>'
)


# Write navbar.html
dir.create("_includes")
writeLines(navbar_html, "_includes/navbar.html")
navbar.path = normalizePath("_includes/navbar.html", mustWork = TRUE)


# :: render html ---------------------------------------------------------


# List files
files <- c("index.Rmd", "about.Rmd", "mlp_validation/mlp_validation.Rmd",
           "ml_models/ml_figures.Rmd")

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
  
  # Conditionally add navbar
  if (!(basename(f) %in% c("index.Rmd", "about.Rmd"))) {
    output_options$includes <- list(before_body = navbar.path)
  }
  
  # Render the file
  rmarkdown::render(f, 
                    output_dir = "docs",
                    output_format = "html_document", 
                    output_options = output_options)
}



# Script to render the Rmarkdown file as a github markdown document
rmarkdown::render("Home-Field-Advantage-in-the-Summer-Olympics.Rmd", 
                  output_format = "github_document", 
                  output_dir = "./", 
                  output_options = list(html_preview = FALSE))
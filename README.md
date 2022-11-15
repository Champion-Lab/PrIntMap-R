# PrIntMap-R
A shiny app for viewing protein coverage in a quantitative manner  

Champion Lab  
University of Notre Dame  
By Simon D. Weaver and Christine M. DeRosa  

## Protein Intensity Mapper

Online version can be found hosted on shiny apps:  

https://championlab.shinyapps.io/printmap-r/

Latest Release can be found at Zenodo: 
[![DOI](https://zenodo.org/badge/497019874.svg)](https://zenodo.org/badge/latestdoi/497019874)


To run locally in R:  

 * Download and install [R](https://www.r-project.org/) and optionally (but encouraged), [RStudio](https://www.rstudio.com/products/rstudio/download/).
 * Either clone this directory, use the latest release from Zenodo (above), or download a zip file containing the contents of this repository by navigating to [https://github.com/Champion-Lab/PrIntMap-R/zipball/master](https://github.com/Champion-Lab/PrIntMap-R/zipball/master). If necessary, unzip this file. The location of this on your local computer is the 'App directory'.
 * Open R (or RStudio if installed).
 * Ensure that the shiny package is loaded by running `install.packages("shiny")` and active `library(shiny)` in the R console.
 * Ensure that all of the packages you need are installed on your computer. The required packages can be found in 'global.R'. You can either install them manually, or run this code in the R console:
 ```
 install.packages(c("stringr", "ggplot2", "seqinr", "plotly", "scales", "dplyr", "tidyr", "shinycssloaders", "markdown", "data.table", "readr", "purrr", "shinyBS"))
 ```  
 * Change your working directory to the path of the cloned or downloaded app directory `setwd('path/to/directory')`
 * In Rstudio, open either 'ui.R' or 'Server.R'.
 * A 'Run App' button should appear in the top right of the '.R' file that you just opened. Click this to start a local version of the app. 
 * Alternatively, you can launch the app directly from the R console: `runApp("path/to/appDirectory")`. If you already set the working directory to this path, you can simply run: `runApp()`.
 * If there are errors, you may not have all the correct packages installed. Open the 'global.R' file and ensure that all the packages used are in your 'packages' tab.

# PrIntMap-R
A shiny app for viewing protein coverage in a quantitative manner  

Simon D. Weaver  
Champion Lab  
University of Notre Dame  

## Protein Intensity Mapper

Beta version can be found hosted on shiny apps:  

https://championlab.shinyapps.io/printmap-r/


To run locally in RStudio:
 * Download and install [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/)
 * Ensure that the shiny package is loaded `install.packages(shiny)` and active `library(shiny)`
 * Change your working directory to the path that contains 'ui.R' and 'Server.R' `setwd('path/to/directory')`
 * Open either 'ui.R' or 'Server.R'
 * A 'Run App' button should appear in the top right of the '.R' file that you just opened. Click this to start a local version of the app
 * If there are errors, you may not have all the correct packages downloaded. Open the 'global.R' file and ensure that all the packages used are in your 'packages' tab. Alternatively you can simply install all of them with `install.packages()`.

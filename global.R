#load packages
library(stringr)
library(ggplot2)
library(seqinr)
library(plotly)
library(scales)
library(dplyr)
library(tidyr)
library(shinycssloaders)
library(markdown)
library(data.table)
library(readr)
library(purrr)

#reactlog
library(reactlog)
options(shiny.reactlog = TRUE)

#source functions and libraries
source("www/libraries.R")
source("www/basic_functions.R")
source("www/two_sample_functions.R")
source("www/annotation_functions.R")
source("www/stacked_plot_functions.R")
source("www/error_handling.R")
source("www/unique_peptide_functions.R")

#allow user to see error messages
options(shiny.sanitize.errors = FALSE)


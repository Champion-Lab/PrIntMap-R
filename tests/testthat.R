library(testthat)
library(shinytest2)
# test_dir(
#   "./tests/testthat",
#   env = shiny::loadSupport(),
#   reporter = c("progress", "fail")
# )
shinytest2::test_app()

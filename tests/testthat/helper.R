
# library(rprojroot)

quickLookup <- function(name) {
  rprojroot::find_testthat_root_file("historical_data", name)
}

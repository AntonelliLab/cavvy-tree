dir_gen <- function(script) {
  # Make the results directory
  flpth <- file.path('results', script)
  if (!dir.exists(flpth)) dir.create(flpth)
  flpth
}

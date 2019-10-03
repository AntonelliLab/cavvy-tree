
# dirs in stages/
outdir_get <- function(nm) {
  flpth <- file.path('stages', nm)
  if (!dir.exists(flpth)) dir.create(flpth)
  flpth
}
indir_get <- function(nm) {
  flpth <- file.path('stages', nm)
  if (!dir.exists(flpth)) stop('No ', flpth, ' found.', call. = FALSE)
  flpth
}

# filepaths in data/
dataflpth_get <- function(nm) {
  flpth <- file.path('data', nm)
  if (!file.exists(flpth)) stop('No ', flpth, ' found.', call. = FALSE)
  flpth
}

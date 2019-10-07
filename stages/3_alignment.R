# Local functions ----
outdir_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Library ----
library(outsider)

# Vars ----
repo <- 'dombennett/om..mafft'
indir <- indir_get('2_cluster')
outdir <- outdir_get('3_alignment')

# Install ----
if (!is_module_installed(repo)) {
  module_install(repo = repo, force = TRUE)
}
mafft <- module_import(fname = 'mafft', repo = repo)

# Align ----
input_files <- list.files(path = indir, pattern = '.fasta')
cids <- sub(pattern = '\\.fasta', replacement = '', x = input_files)
cids <- sub(pattern = '^sequences_', replacement = '', x = cids)
for (i in seq_along(input_files)) {
  cat('... aligning ', crayon::green(cids[[i]]), '\n')
  input_file <- file.path(indir, input_files[[i]])
  output_file <- file.path(outdir, paste0('alignment_', cids[[i]], '.fasta'))
  mafft(arglist = c('--auto', input_file, '>', output_file))
}

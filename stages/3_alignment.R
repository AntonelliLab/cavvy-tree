# Local functions ----
outdir_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Library ----
library(outsider)
verbosity_set(show_program = FALSE, show_docker = FALSE)

# Vars ----
alignment_repo <- 'dombennett/om..pasta'
trimming_repo <- 'dombennett/om..trimal'
indir <- indir_get('2_cluster')
outdir <- outdir_get('3_alignment')

# Install ----
if (!is_module_installed(alignment_repo)) {
  module_install(repo = alignment_repo, force = TRUE)
}
pasta <- module_import(fname = 'pasta', repo = alignment_repo)
if (!is_module_installed(trimming_repo)) {
  module_install(repo = trimming_repo, force = TRUE)
}
trimal <- module_import(fname = 'trimal', repo = trimming_repo)

# Align ----
input_files <- list.files(path = indir, pattern = '.fasta')
cids <- sub(pattern = '\\.fasta', replacement = '', x = input_files)
cids <- sub(pattern = '^sequences_', replacement = '', x = cids)
for (i in seq_along(input_files)) {
  cat('... aligning ', crayon::green(cids[[i]]), '\n')
  input_file <- file.path(indir, input_files[[i]])
  output_file <- file.path(outdir, paste0('alignment_', cids[[i]], '.fasta'))
  pasta(arglist = c('-i', input_file, '-d', 'dna',
                    '--alignment-suffix=alignment.fasta',
                    paste0('--job=', cids[[i]])), outdir = outdir)
}

# Trimming ----
input_files <- list.files(path = outdir, pattern = '\\.alignment\\.fasta')
cids <- sub(pattern = '\\.alignment\\.fasta', replacement = '', x = input_files)
for (i in seq_along(input_files)) {
  cat('... trimming ', crayon::green(cids[[i]]), '\n')
  input_file <- file.path(outdir, input_files[[i]])
  output_file <- paste0('trimmed_alignment_', cids[[i]], '.fasta')
  trimal(arglist = c('-in', input_file, '-out', output_file, '-automated1',
                     '-fasta'), outdir = outdir)
}

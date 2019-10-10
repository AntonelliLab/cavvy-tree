# Local functions ----
outdir_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Functions ----
partitionfile_get <- function(flpth) {
  lns <- readLines(file.path(flpth, 'analysis', 'best_scheme.txt'))
  start_line <- max(which(grepl(pattern = 'raxml', x = lns,
                                ignore.case = TRUE))) + 1
  end_line <- min(which(grepl(pattern = 'mrbayes', x = lns,
                              ignore.case = TRUE))) - 1
  lns <- lns[start_line:end_line]
  lns <- lns[lns != '']
  output_file <- file.path(flpth, 'raxml_partition.txt')
  writeLines(text = lns, con = output_file)
  output_file
}

# Library ----
library(outsider)
verbosity_set(show_program = TRUE, show_docker = FALSE)

# Vars ----
repo <- 'dombennett/om..raxml'
indir <- indir_get('4_supermatrix')
outdir <- outdir_get('5_phylogeny')

# Install ----
if (!is_module_installed(repo)) {
  module_install(repo = repo, force = TRUE)
}
raxml <- module_import(fname = 'raxml', repo = repo)

# Phylogeny ----
supermatrix_files <- list.dirs(path = file.path(indir,
                                                'partitionfinder_results'),
                               recursive = FALSE)
for (i in seq_along(supermatrix_files)) {
  nm <- basename(supermatrix_files[[i]])
  results_dir <- file.path(outdir, nm)
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
  cat('... tree for ', crayon::green(nm), ' supermatrix\n')
  supermatrix_file <- file.path(supermatrix_files[[i]], 'supermatrix.phy')
  partition_file <- partitionfile_get(flpth = supermatrix_files[[i]])
  # fast analysis
  seed_n <- round(runif(n = 1, min = 1, max = 99999))
  raxml(arglist = c('-m', 'GTRGAMMA', '-s', supermatrix_file, '-p', seed_n, '-n',
                    nm, '-T', '4', '-q', partition_file),
        outdir = results_dir)
}

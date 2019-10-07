# Local functions ----
outdir_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Library ----
library(outsider)

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
supermatrix_file <- file.path(indir, 'supermatrix.fasta')
partition_file <- file.path(indir, 'partition.txt')
# fast analysis
seed_n <- round(runif(n = 1, min = 1, max = 99999))
raxml(arglist = c('-m', 'GTRGAMMA', '-s', supermatrix_file, '-p', seed_n, '-n',
                  'supermatrix', '-T', '2', '-q', partition_file),
      outdir = outdir)

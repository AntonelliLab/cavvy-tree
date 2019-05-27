# Local functions ----
source(file.path('tools', 'pipeline.R'))

# Libs ----
library(restez)
library(phylotaR)

# Restez setup ----
restez_path_set('data')

# Run ----
outdir <- dir_gen('1_phylotar')
# for linux
ncbi_dr <- file.path("bin", "ncbi-blast-2.7.1+", "bin")
# Hystricomorpha
txid <- 33550
setup(wd = outdir, txid = txid, ncbi_dr = ncbi_dr, overwrite = TRUE)
run(wd = outdir)
#restart(wd = outdir)

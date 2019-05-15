# Local functions ----
source(file.path('tools', 'pipeline.R'))

# Run ----
outdir <- dir_gen('1_phylotar')
# for linux
ncbi_dr <- file.path("bin", "ncbi-blast-2.7.1+", "bin")
# Hystricomorpha
txid <- 33550
phylotaR::setup(wd = outdir, txid = txid, ncbi_dr = ncbi_dr, overwrite = TRUE)
phylotaR::run(wd = outdir)

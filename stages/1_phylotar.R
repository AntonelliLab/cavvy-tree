# Local functions ----
dataflpth_get <- outdir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Libs ----
library(restez)
library(phylotaR)
outsider::verbosity_set(show_program = FALSE, show_docker = FALSE)

# Restez setup ----
restez_path_set('data')
restez_connect()
if (!restez_ready()) {
  stop('restez not ready.')
}

# Names ----
cavvies_tbl <- read.csv(dataflpth_get('caviomorpha.csv'))
ids <- cavvies_tbl[['Ids']]

# Run ----
outdir <- outdir_get('1_phylotar')
setup(wd = outdir, txid = ids, overwrite = TRUE, outsider = TRUE, ncps = 2,
      btchsz = 200, db_only = FALSE)
run(wd = outdir)
#restart(wd = outdir)

# Clean-up ----
restez_disconnect()

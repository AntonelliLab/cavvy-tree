# Local functions ----
dataflpth_get <- outdir_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Functions ----
duplicates_drop <- function(phylota) {
  # Correct phylotaR error: duplicated clusters for multiple ID search
  seeds <- phylotaR::get_clstr_slot(phylota = phylota, cid = phylota@cids,
                                    slt_nm = 'seed')
  dups <- unique(seeds[duplicated(seeds)])
  to_drop <- NULL
  for (dup in dups) {
    cids <- phylota@cids[seeds == dup]
    nsqs <- phylotaR::get_clstr_slot(phylota = phylota, cid = cids,
                                     slt_nm = 'nsqs')
    to_drop <- c(to_drop, cids[-1 * which.max(nsqs)])
  }
  keep <- phylota@cids[!phylota@cids %in% to_drop]
  res <- drop_clstrs(phylota = phylota, cid = keep)
  res
}

# Libs ----
library(phylotaR)

# Vars ----
min_nspp <- 5
min_mad <- 0.4
indir <- indir_get('1_phylotar')
outdir <- outdir_get('2_cluster')

# Input ----
all_cls <- duplicates_drop(read_phylota(indir))

# Reduce ----
n_taxa <- get_clstr_slot(phylota = all_cls, cid = all_cls@cids, slt_nm = 'ntx')
# drop
keep <- names(n_taxa)[n_taxa > 10]
all_cls <- drop_clstrs(phylota = all_cls, cid = keep)
# reduce to species
species_cls <- drop_by_rank(phylota = all_cls, rnk = 'species', n = 1,
                            choose_by = c("pambgs", "age", "nncltds"),
                            greatest = c(FALSE, FALSE, TRUE))
smmry <- summary(species_cls)
smmry <- smmry[smmry[['MAD']] > min_mad, ]
smmry <- smmry[smmry[['N_taxa']] > min_nspp, ]
slctd <- drop_clstrs(phylota = species_cls, as.character(smmry$ID))
write.csv(smmry, file.path(outdir, 'clusters_summary.csv'))

# Output ----
not_dropped <- dropped <- NULL
for (i in seq_along(slctd@cids)) {
  cid <- slctd@cids[i]
  sids <- slctd@clstrs[[cid]]@sids
  txids <- get_txids(slctd, cid = cid, rnk = 'species')
  scnms <- get_tx_slot(slctd, txids, 'scnm')
  scnms <- gsub('\\s+', '_', scnms)
  to_drop <- grepl(pattern = '([0-9]|\\.)', x = scnms)
  if (sum(!to_drop) < 5) {
    next
  }
  dropped <- c(dropped, scnms[to_drop])
  not_dropped <- c(not_dropped, scnms[!to_drop])
  infile <- file.path(outdir, paste0('sequences_', cid, '.fasta'))
  write_sqs(phylota = slctd, outfile = infile, sid = sids[!to_drop],
            sq_nm = scnms[!to_drop])
}

# Stats ----
dropped <- sort(unique(dropped))
not_dropped <- sort(unique(not_dropped))
species_summary <- data.frame(sp = c(not_dropped, dropped),
                              kept = c(rep(1, length(not_dropped)),
                                       rep(0, length(dropped))))
write.csv(smmry, file.path(outdir, 'species_summary.csv'))
cat(crayon::red(length(not_dropped)), ' species represented!\n', sep = '')


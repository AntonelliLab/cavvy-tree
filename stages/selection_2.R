library(phylotaR)

# INPUT
all_cls <- read_phylota(file.path('results', '1_phylotar'))
ids_of_interest <- read.csv(file = file.path('data', 'caviomorpha.csv'),
                            header = TRUE, stringsAsFactors = FALSE)
genera_ids <- ids_of_interest[ ,3]
genera_nms <- ids_of_interest[ ,2]

# REDUCE TO TAXA OF INTEREST
# IDENTIFY IDS
all_txids <- NULL
for (gnra_id in genera_ids) {
  all_txids <- c(all_txids, phylotaR:::descendants_get(id = gnra_id,
                                                       txdct = all_cls@txdct))
}
# manually add singletons
# "Dinomys" - 108858
# "Trichys" - 73865
# "Lagostomus" - 10154
# "Cuscomys" - 1567517
# "Octomys" - 135583
# "Pipanacoctomys" - 227730
# "Salinoctomys" - 1567519 (no sequences)
# "Spalacopus" - 61880
# "Octodontomys" - 170740
# "Callistomys" - 1567507
# "Carterodon" - 1567513
# "Diplomys" - 1567522
# "Echimys" - 30621
# "Euryzygomatomys" - 43325
# "Hoplomys" - 176501
# "Kannabateomys" - 176503
# "Lonchothrix" - 176505
# "Myocastor" - 10157
# "Olallamys" - 1400528
# "Santamartamys" - 1568970
sngltns <- c(108858, 73865, 10154, 1567517, 135583, 227730, 1567519, 61880,
             170740, 1567507, 1567513, 1567522, 30621, 43325, 176501,
             176503, 176505, 10157, 1400528, 1568970)
all_txids <- c(all_txids, genera_ids)
# all_txids <- c(all_txids, genera_ids, sngltns)
# drop clusters
all_sids_txids <- get_txids(phylota = all_cls, sid = all_cls@sids)
sids_keep <- all_cls@sids[all_sids_txids %in% all_txids]
fltrd_cls <- all_cls
for (i in seq_along(all_cls@cids)) {
  cid <- all_cls@cids[[i]]
  fltrd_cls <- drop_sqs(phylota = fltrd_cls, cid = cid, sid = sids_keep)
}
nsqs <- get_clstr_slot(fltrd_cls, all_cls@cids, slt_nm = 'nsqs')
cids <- all_cls@cids[nsqs > 3]
fltrd_cls <- drop_clstrs(phylota = fltrd_cls, cid = cids)

# REDUCE TO JUST SPECIES
spp_cls <- drop_by_rank(fltrd_cls, rnk = 'species', n = 1,
                        choose_by = c("pambgs", "age", "nncltds"),
                        greatest = c(FALSE, FALSE, TRUE))
# summary
smmry <- summary(spp_cls)
res <- lapply(phylota@cids, get_row)
res <- matrix(unlist(res), nrow = length(phylota@cids), byrow = TRUE)
colnames(res) <- c("ID", "Type", "Seed", "Parent", "N_taxa", 
                   "N_seqs", "Med_sql", "MAD", "Definition", "Feature")
res <- data.frame(res, stringsAsFactors = FALSE)
res[["N_taxa"]] <- as.integer(res[["N_taxa"]])
res[["N_seqs"]] <- as.integer(res[["N_seqs"]])
res[["Med_sql"]] <- as.numeric(res[["Med_sql"]])
res[["MAD"]] <- as.numeric(res[["MAD"]])
res


smmry <- smmry[smmry[['MAD']] > 0.5, ]
smmry <- smmry[smmry[['N_taxa']] > 10, ]
smmry <- smmry[smmry$Type == 'merged', ]
slctd <- drop_cls(spp_cls, as.character(smmry$ID))
write.csv(smmry, file.path('results', 'summary.csv'))

# OUTPUT
# write out top 10 clusters with most taxa
mssng_table <- matrix(NA, nrow = length(slctd@cids),
                      ncol = length(genera_nms))
colnames(mssng_table) <- genera_nms
for (i in seq_along(slctd@cids)) {
  cid <- slctd@cids[i]
  sids <- slctd@cls[[cid]]@sids
  txids <- get_txids(slctd, cid = cid, rnk = 'species')
  scnms <- get_tx_slot(slctd, txids, 'scnm')
  mssng_table[i, ] <- sapply(genera_nms, function(pttrn) any(grepl(pttrn, x = scnms)))
  # filter out non-described species
  pull <- !grepl(pattern = '[^a-zA-Z ]', x = scnms)
  infile <- file.path('phylogeny_files', paste0('sequences_', cid, '.fasta'))
  write_sqs(phylota = slctd, outfile = infile, sid = sids[pull],
            sq_nm = scnms[pull])
}
colnames(mssng_table)[colSums(mssng_table) == 0]

# ALIGN
for (cid in slctd@cids) {
  inpt <- file.path('phylogeny_files', paste0('sequences_', cid, '.fasta'))
  otpt <- file.path('phylogeny_files', paste0('alignment_', cid,'.fasta'))
  system(paste0('mafft --auto ', inpt, ' > ', otpt))
}

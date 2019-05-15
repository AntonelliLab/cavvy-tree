library(phylotaR)

# INPUT
all_cls <- read_phylota('phylotar')

# drop
n_taxa <- get_cl_slot(all_cls, cid = all_cls@cids, slt_nm = 'ntx')
keep <- names(n_taxa)[n_taxa > 10]
all_cls <- drop_cls(all_cls, keep)
table(sapply(all_cls@cls@cls, function(x) x@typ))

# reduce to species
genus_cls <- drop_by_rank(all_cls, rnk = 'genus', n = 2,
                          choose_by = c("pambgs", "age", "nncltds"),
                          greatest = c(FALSE, FALSE, TRUE))
smmry <- summary(genus_cls)
smmry <- smmry[smmry[['MAD']] > 0.5, ]
smmry <- smmry[smmry[['N_taxa']] > 25, ]
smmry <- smmry[smmry$Type == 'merged', ]
smmry <- smmry[order(smmry$N_taxa, decreasing = TRUE)[1:10], ]
slctd <- drop_cls(genus_cls, as.character(smmry$ID))
write.csv(smmry, file.path('results', 'summary.csv'))

# OUTPUT
# write out top 10 clusters with most taxa
for (i in seq_along(slctd@cids)) {
  cid <- slctd@cids[i]
  sids <- slctd@cls[[cid]]@sids
  txids <- get_txids(slctd, cid = cid, rnk = 'genus')
  scnms <- get_tx_slot(slctd, txids, 'scnm')
  n <- sapply(seq_along(scnms), function(x) 
    sum(scnms[x] == scnms[x:length(scnms)]))
  sq_nm <- paste0(scnms, '_', n)
  infile <- file.path('phylogeny_files', paste0('sequences_', cid, '.fasta'))
  write_sqs(phylota = slctd, outfile = infile, sid = sids,
            sq_nm = sq_nm)
}

# ALIGN
for (cid in slctd@cids) {
  inpt <- file.path('phylogeny_files', paste0('sequences_', cid, '.fasta'))
  otpt <- file.path('phylogeny_files', paste0('alignment_', cid,'.fasta'))
  system(paste0('mafft --auto ', inpt, ' > ', otpt))
}


# drop irrevelant sequences
# octodontidae <- 10158
# chinchillidae <- 10150
# ids_of_interest <- c(phylotaR:::getADs(id = octodontidae, txdct = all_cls@txdct),
#                      phylotaR:::getADs(id = chinchillidae, txdct = all_cls@txdct))
# all_sids_txids <- get_txids(phylota = all_cls, sid = all_cls@sids)
# sids_keep <- all_cls@sids[all_sids_txids %in% ids_of_interest]
# fltrd_cls <- all_cls
# for (i in seq_along(all_cls@cids)) {
#   cid <- all_cls@cids[[i]]
#   fltrd_cls <- drop_sqs(phylota = fltrd_cls, cid = cid, sid = sids_keep)
# }
# # drop empty clusters
# nsqs <- get_cl_slot(fltrd_cls, all_cls@cids, slt_nm = 'nsqs')
# cids <- all_cls@cids[nsqs > 3]
# fltrd_cls <- drop_cls(phylota = fltrd_cls, cid = cids)
# # summary
# smmry <- summary(fltrd_cls)
# smmry <- smmry[smmry[['MAD']] > 0.5, ]
# smmry <- smmry[smmry[['N_taxa']] > 10, ]


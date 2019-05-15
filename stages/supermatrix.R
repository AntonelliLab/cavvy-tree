# Combine alignments into a supermatrix, run RAxML

# LIBS
source(file.path('tools', 'supermatrix.R'))

# INPUT
alfls <- list.files('phylogeny_files', pattern='alignment_[0-9]+.fasta')
als <- vector('list', length = length(alfls))
for (i in seq_along(alfls)) {
  als[[i]] <- readSqs(file.path('phylogeny_files', alfls[i]))
}

# DROP OVERHANGING EDGES
(sapply(als, function(x) nchar(x[[1]])))
als <- drpOvrhngs(als, ctff=0.75)
(sapply(als, function(x) nchar(x[[1]])))
n_taxa <- sapply(als, length)

# GEN PARTITION TEXT
lngths <- sapply(als, function(x) nchar(x[[1]]))
partition(lngths, fl=file.path('phylogeny_files', 'partition.txt'))

# GEN SUPERMARTIX
fllrs <- sapply(lngths, function(x) paste0(rep('-', x),
                                           collapse = ''))
all_nms <- unique(unlist(sapply(als, names)))
all_nms <- sort(all_nms)
pull <- !grepl('\\ssp\\.', all_nms)
all_nms <- all_nms[pull]
sprmtrx <- vector('list', length = length(all_nms))
names(sprmtrx) <- all_nms
for (nm in all_nms) {
  al <- ''
  for (i in seq_along(als)) {
    tmp <- als[[i]][[nm]]
    tmp <- ifelse(is.null(tmp), fllrs[[i]], tmp)
    al <- paste0(al, tmp)
  }
  sprmtrx[[nm]] <- al
}

# DROP TAXA WITH TOO MANY GAPS
ngaps <- sapply(gregexpr('-', sprmtrx), length)
pull <- ngaps < nchar(sprmtrx[[1]])
sprmtrx <- sprmtrx[pull]

# CHECK AND WRITE OUT
all(sapply(sprmtrx, nchar) == nchar(sprmtrx[[1]]))
names(sprmtrx) <- gsub('\\s', '_', names(sprmtrx))
writeSqs(sprmtrx, fl = file.path('phylogeny_files', 'supermatrix.fasta'))

# RAxML
# Warning: partition.txt may need minor modification depending on gene type
inpt <- file.path('phylogeny_files', 'supermatrix.fasta')
prttnfl <- file.path('phylogeny_files', 'partition.txt')
system(paste0('raxmlHPC -f a -m GTRGAMMA -T 2 -# 10 -p ',
              sample(0:10000000, 1), ' -x ', sample(0:10000000, 1),
              ' -n test2 -s ', inpt, ' -q ', prttnfl))
# consensus
system('raxmlHPC -m GTRCAT -J MR -z RAxML_bootstrap.palms -n palms_con')

# CLEAN-UP
file.rename('RAxML_bestTree.palms', file.path(wd, 'best_tree.tre'))
file.rename('RAxML_bootstrap.palms', file.path(wd, 'bootstraps.tre'))
file.rename('RAxML_MajorityRuleConsensusTree.palms_con',
            file.path(wd, 'consensus.tre'))
file.remove(list.files(pattern = 'RAxML.*\\.palms'))
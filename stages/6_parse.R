# Local functions ----
dataflpth_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Library ----
library(ape)

# Vars ----
indir <- indir_get('5_phylogeny')
input_file <- file.path(indir, 'low', 'RAxML_bestTree.low')

# Read in ----
tree <- read.tree(file = input_file)
# root tree
chinchilloidea_pattern <- '(Hystrix|Atherurus|Trichys|Lagidium|Chinchilla|Dinomys|Lagostomus)'
outgroup_tips <- tree$tip.label[grepl(pattern = chinchilloidea_pattern,
                                      x = tree$tip.label, ignore.case = TRUE)]
tree <- root(tree, outgroup = outgroup_tips, resolve.root = TRUE)

# Stats ----
cavvies_tbl <- read.csv(dataflpth_get('caviomorpha.csv'),
                        stringsAsFactors = FALSE)
genus_presence <- vapply(X = cavvies_tbl[['Genus']], FUN = function(x) {
  any(grepl(pattern = x, x = tree$tip.label, ignore.case = TRUE))
  }, FUN.VALUE = logical(1))
sum(genus_presence) == length(genus_presence)
missing_genera <- names(genus_presence)[!genus_presence]

# Plot ----
msg <- paste0(length(tree$tip.label), ' Nspp. | ')
png(filename = 'tree.png', units = 'px', width = 480*1.75, height = 480*1.75)
if (length(missing_genera) > 0) {
  msg <- paste0(msg, 'Missing genera: ', paste0(missing_genera,
                                                collapse = ', '))
} else {
  msg <- paste0(msg, 'No missing genera')
}
par(mar = c(1, 0, 0, 0))
plot(x = tree, type = 'fan')
mtext(text = msg, side = 1, line = -1, cex = 0.75)
dev.off()

# Write out ----
readme_txt <- '
# Cavvy-tree

Results for *cavvy-tree* project: "construction of phylogenetic supermatrix"
for select Caviomorph families.

Three supermatrices are made available for different strictnesses of
supermatrix-completeness. "low" provides the matrix with maximum species
coverage, but lots of gaps. "high" provides minimum number of gaps but has
fewest number of species.

Matrices are provided in FASTA and PHYLIP format.

phylotaR results are summarised in "clusters_summary" and "species_summary".

For each matrix type a RAxML guide tree and best partition schemes according to
PartitionFinder2 are provided.

For more information visit: https://github.com/AntonelliLab/cavvy-tree

## Folder structure

```
results/
- clusters_summary.csv
- species_summary.csv
- - low/
- - - RAxML_bestTree
- - - supermatrix.fasta
- - - supermatrix.phy
- - - best_scheme.txt
- - mid/
- - - RAxML_bestTree
- - - supermatrix.fasta
- - - supermatrix.phy
- - - best_scheme.txt
- - high/
- - - RAxML_bestTree
- - - supermatrix.fasta
- - - supermatrix.phy
- - - best_scheme.txt
- README.txt
```
'
readme_txt <- paste0(readme_txt, '\n\n**Date and time**: ', '\n', Sys.time())
if (!dir.exists('results')) {
  dir.create('results')
}
files_to_copy <- file.path('stages', '2_cluster', c('clusters_summary.csv',
                                                    'species_summary.csv'))
file.copy(from = files_to_copy, to = 'results', overwrite = TRUE)
analysis_dirs <- file.path('results', c('low', 'mid', 'high'))
for (analysis_dir in analysis_dirs) {
  if (!dir.exists(analysis_dir)) {
    dir.create(analysis_dir)
  }
  nm <- basename(analysis_dir)
  files_to_copy <- c(file.path(indir, nm, paste0('RAxML_bestTree.', nm)),
                     file.path('stages', '4_supermatrix',
                               'partitionfinder_results', nm,
                               c('supermatrix.fasta', 'supermatrix.phy')),
                     file.path('stages', '4_supermatrix',
                               'partitionfinder_results', nm, 'analysis',
                               'best_scheme.txt'))
  file.copy(from = files_to_copy, to = analysis_dir, overwrite = TRUE)
}
cat(readme_txt, file = file.path('results', 'README.md'))

# Compress ----
utils::zip(zipfile = 'results.zip', files = 'results')

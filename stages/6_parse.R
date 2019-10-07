# Local functions ----
indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Library ----
library(ape)

# Vars ----
indir <- indir_get('5_phylogeny')
input_file <- file.path(indir, 'RAxML_bestTree.supermatrix')

# Read in ----
tree <- read.tree(file = input_file)
# drop repeated tips
pattern <- '_[0-9]$'
to_drop <- tree$tip.label[duplicated(sub(pattern = pattern,
                                         x = tree$tip.label,
                                         replacement = ''))]
tree <- drop.tip(phy = tree, tip = to_drop)
tree$tip.label <- sub(pattern = pattern, x = tree$tip.label, replacement = '')
# root tree
outgroup_tips <- tree$tip.label[!grepl(pattern = '(Hystrix|Atherurus|Trichys)',
                                      x = tree$tip.label, ignore.case = TRUE)]
tree <- root(tree, outgroup = outgroup_tips, resolve.root = TRUE)

# Plot ----
png(filename = 'tree.png', units = 'px', width = 480*1.75, height = 480*1.75)
plot(x = tree, no.margin = TRUE, type = 'fan')
dev.off()

# Write out ----
if (!dir.exists('results')) {
  dir.create('results')
}
write.tree(phy = tree, file = file.path('results', 'guide_tree.tre'))
sprmtrx_flpth <- file.path('stages', '4_supermatrix')
fls <- file.path(sprmtrx_flpth, list.files(sprmtrx_flpth))
file.copy(from = fls, to = 'results', overwrite = TRUE)

# Compress ----
utils::zip(zipfile = 'results.zip', files = 'results')

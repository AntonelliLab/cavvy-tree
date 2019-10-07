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
to_drop <- tree$tip.label[!duplicated(sub(pattern = '_[0-9]$',
                                          x = tree$tip.label,
                                          replacement = ''))]
tree <- drop.tip(phy = tree, tip = to_drop)
tree$tip.label <- sub(pattern = '_[0-9]$', x = tree$tip.label, replacement = '')
# # root tree
# outgroup_tips <- tree$tip.label[grepl(pattern = outgroup_pattern,
#                                       x = tree$tip.label, ignore.case = TRUE)]
# tree <- root(unroot(tree), outgroup = outgroup_tips)
png(filename = 'tree.png', units = 'px', width = 480*1.75, height = 480*1.75)
plot(x = tree, no.margin = TRUE, type = 'fan')
dev.off()


# Local functions ----
outdir_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Library ----
library(gaius)

# Vars ----
indir <- indir_get('3_alignment')
outdir <- outdir_get('4_supermatrix')

# Input ---
alignment_files <- file.path(indir, list.files(path = indir,
                                               pattern = '.fasta'))
alignment_list <- alignment_read(flpths = alignment_files)

# Supermatrix ----
# min_ngenes - minimum number of genes in matrix
# min_ntips - minimum number of tips in matrix
# min_nbps - minimum number of base paris in a gene
# column_cutoff - proportion of non-gaps per column
# tip_cutoff - proportion of non-gaps per tip
supermatrix <- supermatrices_get(alignment_list = alignment_list,
                                 min_ngenes = 2, min_ntips = 5, min_nbps = 250,
                                 column_cutoff = 0.1, tip_cutoff = 0.1)
supermatrix <- supermatrix[[1]]
# number of bps per cluster
nbps <- attr(supermatrix, 'nbps')

# Write out ----
partition_file(nbp = nbps, flpth = file.path(outdir, "partition.txt"))
sequences_write(x = supermatrix, flpth = file.path(outdir, 'supermatrix.fasta'))

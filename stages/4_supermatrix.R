# Local functions ----
outdir_get <- indir_get <- NULL
source(file.path('tools', 'pipeline.R'))

# Functions ----
partition_write <- function(nbps, flpth, wobble = TRUE) {
  # write partition file in PartitionFinder format
  gene <- strt <- 1
  data_block <- NULL
  for (i in seq_along(nbps)) {
    end <- nbps[[i]] + strt - 1
    if (wobble) {
      data_block <- c(data_block, paste0('gene_', gene, "_pos1 = ", 
                                         strt, "-", end, "\\3;"))
      data_block <- c(data_block, paste0('gene_', gene, "_pos2 = ", 
                                         strt + 1, "-", end, "\\3;"))
      data_block <- c(data_block, paste0('gene_', gene, "_pos3 = ", 
                                         strt + 2, "-", end, "\\3;"))
    } else {
      data_block <- c(data_block, paste0('gene_', gene, " = ", 
                                         strt, "-", end, ";"))
    }
    strt <- end + 1
    gene <- gene + 1
  }
  pft <- readLines(con = file.path('data', 'partition_finder_template.cfg'))
  i <- which(pft == '%data_blocks%')
  new_pft <- c(pft[seq_len(i - 1)], data_block,
               pft[seq(from = i + 1, to = length(pft), by = 1)])
  writeLines(text = new_pft, flpth)
}
partition_finder_setup <- function(supermatrix, outdir, nm, wobble) {
  # setup a folder for PartitionFinder
  supermatrix <- supermatrix[[1]]
  nbps <- attr(supermatrix, 'nbps')
  mtrx_dir <- file.path(outdir, nm)
  if (!dir.exists(mtrx_dir)) {
    dir.create(mtrx_dir)
  }
  sequences_write(x = supermatrix, flpth = file.path(mtrx_dir,
                                                     'supermatrix.fasta'))
  partition_write(nbps = nbps, flpth = file.path(mtrx_dir,
                                                 "partition_finder.cfg"),
                  wobble = wobble)
}

# Library ----
library(gaius)
library(outsider)
verbosity_set(show_program = TRUE, show_docker = FALSE)

# Vars ----
pf_repo <- 'dombennett/om..partitionfinder'
readal_repo <- 'dombennett/om..trimal'
indir <- indir_get('3_alignment')
outdir <- outdir_get('4_supermatrix')

# Install ----
if (!is_module_installed(pf_repo)) {
  module_install(repo = pf_repo, force = TRUE)
}
prtfnd <- module_import(fname = 'PartitionFinder', repo = pf_repo)
if (!is_module_installed(readal_repo)) {
  module_install(repo = readal_repo, force = TRUE)
}
readal <- module_import(fname = 'readal', repo = readal_repo)

# Input ---
alignment_files <- file.path(indir, list.files(path = indir,
                                               pattern = 'trimmed_alignment'))
alignment_list <- alignment_read(flpths = alignment_files)

# Supermatrix ----
# vars:
# min_ngenes - minimum number of genes in matrix
# min_ntips - minimum number of tips in matrix
# min_nbps - minimum number of base paris in a gene
# column_cutoff - proportion of non-gaps per column
# tip_cutoff - proportion of non-gaps per tip
# three supermatrices:
# low quality (80% gaps, all data)
supermatrix <- supermatrices_get(alignment_list = alignment_list,
                                 min_ngenes = 2, min_ntips = 5, min_nbps = 200,
                                 column_cutoff = 0.0, tip_cutoff = 0.0)
partition_finder_setup(supermatrix = supermatrix, outdir = outdir, nm = 'low',
                       wobble = TRUE)
# mid quality (50% gaps, columns dropped)
supermatrix <- supermatrices_get(alignment_list = alignment_list,
                                 min_ngenes = 2, min_ntips = 5, min_nbps = 200,
                                 column_cutoff = 0.1, tip_cutoff = 0.3)
# Note: no wobble base pair if columns are dropped
partition_finder_setup(supermatrix = supermatrix, outdir = outdir, nm = 'mid',
                       wobble = FALSE)
# high quality (20% gaps, columns dropped)
supermatrix <- supermatrices_get(alignment_list = alignment_list,
                                 min_ngenes = 2, min_ntips = 5, min_nbps = 200,
                                 column_cutoff = 0.2, tip_cutoff = 0.7)
partition_finder_setup(supermatrix = supermatrix, outdir = outdir, nm = 'high',
                       wobble = FALSE)

# PartitionFinder ----
matrix_dirs <- file.path(outdir, c('high', 'mid', 'low'))
pfdir <- file.path(outdir, 'partitionfinder_results')
if (!dir.exists(pfdir)) {
  dir.create(pfdir)
}
for (i in seq_along(matrix_dirs)) {
  cat('... finding partitions for ', crayon::green(basename(matrix_dirs[[i]])),
      ' supermatrix\n')
  fasta_file <- file.path(matrix_dirs[[i]], 'supermatrix.fasta')
  readal(arglist = c('-in', fasta_file, '-out', 'supermatrix.phy',
                     '-phylip'), outdir = matrix_dirs[[i]])
  prtfnd(arglist = matrix_dirs[[i]], outdir = pfdir)
}

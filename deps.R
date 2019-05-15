# To successfully run this pipeline, the following packages/modules are required

# R packages ----
devtools::install_github('ropensci/phylotar')
# if fails.... try
# devtools::install_github("hannesmuehleisen/MonetDBLite-R")
# devtools::install_github('ropensci/restez')
devtools::install_github('antonellilab/gaius')
devtools::install_github('antonellilab/outsider.base')
devtools::install_github('antonellilab/outsider')

# outsider mods ----
outsider::module_install(repo = 'dombennett/om..raxml')
outsider::module_install(repo = 'dombennett/om..mafft')

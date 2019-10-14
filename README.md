# Caviomorpha phylogenetic tree of extinct and extant taxa

> :hamster: Building a tree of some small fluffy animals!

Construction supermatrix of biological sequences for phylogenetic analysis for
["Caviomorpha"](https://en.wikipedia.org/wiki/Caviomorpha), families:
Hystricidae, Chinchillidae, Dinomyidae, Abrocomidae, Octodontidae, Ctenomyidae
and Echimyidae.

## Process

1. Download GenBank (`0_restez.R`)
2. Identify clusters (`1_phylotar.R`)
3. Choose clusters (`2_cluster.R`)
4. Align sequences (`3_alignment.R`)
5. Construct supermatrix (`4_supermatrix.R`)
6. Construct ML phylogeny (`5_phylogeny.R`)
7. Parse results for export (`6_parse.R`)

All stage scripts can be found in `stages/`.

## Key packages

* [`restez`](https://github.com/ropensci/restez)
* [`phylotaR`](https://github.com/ropensci/phylotar)
* [`outsider`](https://github.com/AntonelliLab/outsider)
* [`gaius`](https://github.com/AntonelliLab/gaius)

(All packages are part of the
[`supersmartR`](https://github.com/AntonelliLab/supersmartR) project.)

## Key programs

* [pasta](https://github.com/smirarab/pasta)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [PartitionFinder](http://www.robertlanfear.com/partitionfinder/)
* [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html)
* [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK537770/)

## Methods

Detailed method description can be found in this
[Google Document](https://docs.google.com/document/d/1Q-5b3eSMKInqiy3cgw1JfG54VnTLr_WSFtyMBf1UXz8/edit?usp=sharing).

## Results

### Download

[Download `results.zip`](https://github.com/AntonelliLab/cavvy-tree/raw/master/results.zip)

(Read the `README.md` file containted in the zipped folder to learn about the
file content.)

### Preview

![](https://raw.githubusercontent.com/AntonelliLab/cavvy-tree/master/tree.png)

### Target taxon
![](https://upload.wikimedia.org/wikipedia/commons/f/fc/Two_Adult_Guinea_Pigs_%28cropped%29.jpg)

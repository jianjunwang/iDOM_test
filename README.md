
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iDOM

Ecological networks of dissolved organic matter and microorganisms under
global change.

## About

The package includes functions for analyzing the thermal responses, assembly processes, transformation and microbial associations of DOM.

Package maintainer: Jianjun Wang(<jjwang@niglas.ac.cn>)

Developers: Ang Hu(<anghu@hunau.edu.cn>), Fanfan
Meng(<mengfanfan19@mails.ucas.ac.cn>)

## Installation

You can install the development version of **iDOM** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jianjunwang/iDOM")
```

Installing packages in R from zip source.You may have downloaded
**iDOM** in ***zip*** or ***tar.gz*** format. In order to install the
package from a local zip file you just need to call the install.packages
function with arguments repos = NULL and type = “source”.

``` r
install.packages("file_path/package_file_name.extension",repos = NULL, type = "source")
```

## Key functions in iDOM package

[DOM.H2](https://doi.org/10.1101/2021.08.12.456177): Calculates the
network-level specialization of all interacting trophic levels in
DOM-microbe bipartite networks, including full, negative and positive
networks.

**iCTR**: The function of iCTR is relevant to Hu et al. (under review). The manuscript will be available as preprint soon. The code will be updated upon acceptance.

## References

Wang J, Pan F, Soininen J, Heino J, Shen J. 2016. Nutrient enrichment
modifies temperature-biodiversity relationships in large-scale field
experiments. ***Nature Communications*** 7:1-9.
<https://doi.org/10.1038/ncomms13960>

Hu A, Choi M, Tanentzap AJ, Liu J, Jang K-S, Lennon JT, Liu Y, Soininen
J, Lu X, Zhang Y, Shen J, Wang J. 2021. Quantifying microbial
associations of dissolved organic matter under global change.
***bioRxiv***. <https://doi.org/10.1101/2021.08.12.456177>

Hu A, Choi M, Tanentzap AJ, Liu J, Jang K-S, Lennon JT, Liu Y, Soininen J, Lu X, Zhang Y, Shen J, Wang J. 2022. Ecological networks of dissolved organic matter and microorganisms under global change. ***Nature communications***. <https://www.nature.com/articles/s41467-022-31251-1>

Hu A, Jang K-S, Meng F, Stegen J, Tanentzap AJ, Choi M, Lennon JT, Soininen J, and Wang J. 2022. Microbial and Environmental Processes Shape the Link between Organic Matter Functional Traits and Composition. ***Environmental Science & Technology***. <https://pubs.acs.org/doi/abs/10.1021/acs.est.2c01432>

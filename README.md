# CHETAH: a selective, hierarchical cell type identification method for single-cell RNA sequencing
__CHETAH is a package for cell type identification of single-cell RNA-sequencing (scRNA-seq) data.__
Cell types are assigned by correlating the input data to a reference in a hierarchical manner. CHETAH is built to work with scRNA-seq references, but will also work with RNA-seq or micro-array reference datasets.

CHETAH is submitted to Bioconductor, and will soon be available by:
```{r echo=TRUE, eval=FALSE}
## Install BiocManager is neccesary
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install('CHETAH')

# Load the package
library(CHETAH)
```


The development version can be downloaded from github.
Note that install_github does not always install all dependencies,
so please check if all dependencies are available for you.
```{r echo=TRUE, eval=FALSE}
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jdekanter/CHETAH")

## Check if all dependencies are installed
dep <- c('bioDist', 'ggplot2', 'gplots', 'cowplot',
         'dendextend', 'corrplot', 'reshape2', 'plotly', 'grDevices')
pkg_avail <- suppressMessages(sapply(dep, function (pkg) pkg %in% installed.packages()[, "Package"]))

# --- Install dependencies, if neccesary
if(length(dep[!pkg_avail]) > 0) {
  if (!require("BiocManager")) {
    install.packages("BiocManager")
  }
  BiocManager::install(dep[!pkg_avail])
}
# Load the package
library(CHETAH)
```

To get to know the basics of the CHETAH pacakge, please look at the vignette;
```{r echo=TRUE, eval=FALSE}
vignette("CHETAH_introduction")
```

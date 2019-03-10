# CHETAH: a selective, hierarchical cell type identification method for single-cell RNA sequencing
__CHETAH is a package for cell type identification of single-cell RNA-sequencing (scRNA-seq) data.__
Cell types are assigned by correlating the input data to a reference in a hierarchical manner. CHETAH is built to work with scRNA-seq references, but will also work (with limited capabilities) with RNA-seq or micro-array reference datasets.


A pre-print of the article describing CHETAH can be found at: https://www.biorxiv.org/content/10.1101/558908v1

> NOTE: CHETAH is submitted to Bioconductor.
> CHETAH is now build on top of the `SingleCellExperiment` package,
> which requires some adjustment of running CHETAH compared to the previous versions.

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

At a glance: to run chetah on an input count matrix `input_counts` with t-SNE coordinates in `input_tsne`, and a reference count matrix `ref_counts` with celltypes vector `ref_ct`, run:  

```{r glance, echo=TRUE, eval=FALSE}
## Make SingleCellExperiments
reference <- SingleCellExperiment(assays = list(counts = ref_counts),
                                     colData = DataFrame(celltypes = ref_ct))

input <- SingleCellExperiment(assays = list(counts = input_counts),
                              reducedDims = SimpleList(TSNE = input_tsne))

## Run CHETAH
input <- CHETAHclassifier(input = input, ref_cells = reference)

## Plot the classification
PlotCHETAH(input)

## Extract celltypes:
celltypes <- input$celltype_CHETAH
```

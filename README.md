# CHETAH: a selective, hierarchical cell type identification method for single-cell RNA sequencing
__CHETAH is an R package for cell type identification of single-cell RNA-sequencing (scRNA-seq) data.__
Cell types are assigned by correlating the input data to a reference in a hierarchical manner. CHETAH is built to work with scRNA-seq references, but will also work (with limited capabilities) with RNA-seq or micro-array reference datasets.


The article describing CHETAH can be found at: [Nucleic Acids Research](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz543/5521789?searchresult=1).  

> CHETAH is now part of [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/CHETAH.html).  

CHETAH can be installed by running:
```{r echo=TRUE, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CHETAH")
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

#' A SingleCellExperiment on which CHEATHclassifier is run using the \code{headneck_ref}
#' It holds subset of the Melanoma data, from Tirosh et al. (2016), Science (\code{\link{data_mel}})
#'
#' @docType data
#' @format This is a SingleCellExperiment
#' @references Tirosh et al. (2016) Science 6282:189-196
#'
#' @source for the original data: \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056}{GEO}
"input_mel"

#' A SingleCellExperiment with celltypes in the "celltypes" colData.
#' A subset of the Head-Neck data from Puram et al. (2017) Cancer Cell.
#'
#' @docType data
#' @format A list of expression matrices. Each object is named as the cell type of the cells in that matrix.
#' Each matrix has the cell (names) in the colums and the genes in the rows.
#' @references Puram et al. (2017) Cancer Cell 171:1611-1624
#'
#' @source for the original data: \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322}{GEO}
"headneck_ref"

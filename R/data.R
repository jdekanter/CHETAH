#' The result of running t-SNE on the \code{\link{data_mel}} (subset of Melanoma data).
#' This is matrix of two columns (the two dimensions), with the cells (names) in the rows
#'
#' @docType data
#' @format A matirx with the t-SNE coordinates of the melanoma data in 2 columns (for the two dimensions)
#' @references Tirosh et al. (2016) Science 6282:189-196
#'
#' @source for the original data: \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056}{GEO}
"tsne_mel"

#' A subset of the Melanoma data, from Tirosh et al. (2016), Science (\code{\link{data_mel}})
#'
#' @docType data
#' @format This is a sparse Matrix, with the cells in the columns and the genes in the rows
#' @references Tirosh et al. (2016) Science 6282:189-196
#'
#' @source for the original data: \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056}{GEO}
"data_mel"

#' The CHETAH output of the subset of the Melanoma data (\code{\link{data_mel}})
#' , run with the Head-Neck reference (\code{\link{reference_hn}})
#'
#' @docType data
#' @format see \code{\link{CHETAHclassifier}}
#' @references Tirosh et al. (2016) Science 6282:189-196
#'
#' @source for the original data: \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056}{GEO}
"chetah"

#' A subset of the Head-Neck data from Puram et al. (2017) Cancer Cell.
#'
#' @docType data
#' @format A list of expression matrices. Each object is named as the cell type of the cells in that matrix.
#' Each matrix has the cell (names) in the colums and the genes in the rows.
#' @references Puram et al. (2017) Cancer Cell 171:1611-1624
#'
#' @source for the original data: \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322}{GEO}
"reference_hn"

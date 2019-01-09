## -------------------------------------------------------------------------------------
## ----------------------------------- CHETAH ------------------------------------------
## -- CHaracterization of cEll Types Aided by Hierarchical clustering ------------------
## -------------------------------------------------------------------------------------

#' Identification of cell types aided by hierarchical clustering
#'
#' @description
#' CHETAH classifies an input dataset by comparing it to
#' a reference dataset in a stepwise, top-to-bottom fashion.
#' See 'details' for a full explanation.
#' \emph{NOTE: We recommend to use all the default parameters}
#'
#' @param input \strong{required}: an expression matrix (matrix, data.frame or Matrix)
#' of the cells of interest, with cells in the columns and genes in the rows
#' @param ref_cells \strong{required}: either \cr (1) a list of expression matrices,
#' with one matrix for each reference cell type. IMPORTANT:
#' \emph{the list has to be named!} (with the corresponding cell types) \cr
#' (2) one expression matrix: in this case provide the cell type label
#' in the \code{ref_type} variable.
#' @param ref_types \strong{required}: when ref_cells is óne matrix,
#' ref_types must be a vector of the cell types.
#' The names must correspond to the column names of \code{ref_cells}, in the same order!.
#' @param ref_profiles \emph{optional} In case of bulk-RNA seq or micro-arrays,
#' an expression matrix with one (average) reference expression profile per cell type.
#' @param thresh the initial confidence threshold (can be changed later
#' by \code{\link{Classify}})
#' @param gs_method method for gene selection. In every node of the tree: \cr
#' "fc" = quick method: either a fixed amount (\code{n_genes})
#' of genes is selected with the highest fold-change (default), \cr
#' or genes are selected that have a fold-change higher than \code{fc_thresh}
#' (the latter is used when \code{fix_ngenes = FALSE}) . \cr
#' "wilcox": genes are selected based on fold-change (\code{fc_thresh}),
#' percentage of expression (\code{pc_thresh}) and p-values (\code{p_thresh}),
#' p-values are found by the wilcox test
#' @param cor_method the correlation measure: one of:
#' "spearman" (default), "kendall", "pearson", "cosine"
#' @param clust_method the method used for clustering the reference profiles
#' one of the methods from \code{\link[stats]{hclust}}
#' @param clust_dist a distance measure, default: \code{\link[bioDist]{spearman.dist}}
#' @param n_genes The amount of genes used in every step. Only used if
#' \code{fix_ngenes = TRUE}
#' @param pc_thresh when: \emph{gs_method = "wilcox"}, only genes are selected
#' for which more than a \code{pc_tresh} fraction of a reference group of cells
#' express that gene
#' @param p_thresh when: \emph{gs_method = "wilcox" }, only genes are selected
#' that have a p-value < \code{p_thresh}
#' @param fc_thresh when: \emph{gs_method = "wilcox" or gs_method = "fc"
#' AND fix_ngenes = FALSE},
#' only genes are selected that have a log2 fld-change > \code{fc_thresh}
#' between two reference groups. \cr
#' \strong{if this mode is selected, the reference must be in the log2 space}.
#' @param subsample to prevent reference types with a lot of cells to influence
#' the gene selection, subsample types with more that \code{subsample} cells
#' @param fix_ngenes when: \emph{gs_method = "fc"} use a fixed amount of genes
#' for all correlations. when: \emph{gs_method = "wilcox"}
#' use a maximum of genes per step.
#' When \code{fix_ngenes = FALSE & gs_methode = "fc"} \code{fc_thresh}
#' is used to define the fold-change cut-off for gene selection.
#' @param plot.tree Plot the classification tree.
#' @param only_pos \emph{not recommended}: only use genes for a reference type
#' that are higher expressed in that type, than the others in that node.
#' @param print_steps whether the number of genes (postive and negative)
#' per step per ref_cell_type should be printed
#'
#' @return
#' The chetah object:
#' A list containing the following objects:
#' NOTE: all objects that are lists themselfs contain one object
#' per node in the classification tree
#' \itemize{
#'   \item \strong{classification} a named vector: the classified types
#'   with the corresponding names of the input cells
#'   \item \strong{prof_scores} A list with the profile scores
#'   \item \strong{conf_scores} A list with the confidence scores
#'   \item \strong{nodetypes} A list with the cell types under each node
#'   \item \strong{correlations} A list with the correlations of the
#'   input cells to the reference profiles
#'   \item \strong{tree} the hclust object of the classification tree
#'   \item \strong{nodecoor} the coordinates of the nodes of the classification tree
#'   \item \strong{genes} A list per node, containing a list per
#'   reference type with the genes used for the profile scores of that type
#'   \item \strong{parameters} The parameters used
#' }
#' @details
#' CHETAH will hierarchically cluster reference data
#' to produce a classification tree (ct).
#' The Classification will In each node of the ct, CHETAH will
#' divide the classification to the two branches based on gene selections,
#' correlations and calculation of profile and confidence scores.
#' Some input cells will reach the leaf nodes of the ct (the pre-defined cell types),
#' these classifications are called \strong{final types}
#' For other cells, assignment will stop in a node, if no enough evidence is available
#' to assign that cell to one of the branches of that node. These classifications
#' are called \strong{intermediate types}.  \cr
#' @export
#' @importFrom bioDist spearman.dist
#' @importFrom stats as.dendrogram cutree ecdf hclust p.adjust wilcox.test
#' @examples
#' ## Melanoma data from Tirosh et al. (2016) Science
#' input = data_mel
#' ## Head-Neck data from Puram et al. (2017) Cancer Cell
#' reference = reference_hn
#' ls.str(reference) ## a list of expression matrices, group per cell type
#' chetah <- CHETAHclassifier(input = input, ref_cells = reference)
CHETAHclassifier <- function (input,
                              ref_cells =       NULL,
                              ref_types =       NULL,
                              ref_profiles =    NULL,
                              thresh =          0.1,
                              gs_method =       "fc",
                              cor_method =      "spearman",
                              clust_method =    "average",
                              clust_dist =      bioDist::spearman.dist,
                              n_genes =         200,
                              pc_thresh =       0.2,
                              p_thresh =        0.05,
                              fc_thresh =       1.5,
                              subsample =       FALSE,
                              fix_ngenes =      TRUE,
                              plot.tree =       FALSE,
                              only_pos =        FALSE,
                              print_steps =     FALSE) {

    stopifnot(is.logical(plot.tree),
            is.logical(subsample),
            is.logical(only_pos),
            is.numeric(n_genes),
            is.numeric(pc_thresh),
            is.numeric(p_thresh),
            is.numeric(fc_thresh),
            is.numeric(thresh),
            cor_method %in% c("pearson", "kendall", "spearman", "cosine"),
            gs_method %in% c("wilcox","fc"),
            is(clust_dist(matrix(seq_len(4), nrow = 2)), "dist"),
            !(is.null(ref_cells)) | !(is.null(ref_profiles)),
            if (!is.null(ref_types)) is.character(ref_types) | is.vector(ref_types) else TRUE)

    if (!is.list(ref_cells) & is.null(ref_types)) {
        stop("Please provide a vector of the reference cell types in 'ref_types'")
    }
    if (!is.list(ref_cells) & !is.null(ref_types)) {
        if (!all.equal(colnames(ref_cells), names(ref_types))) {
            stop("the names of 'ref_cells' must be excatly the same as the colnames of 'ref_cells'")
        }
    }
    if (!is.null(ref_profiles)) if (sum(is.na(ref_profiles)) != 0) stop("'ref_profiles' cannot hold NA values'")
    if (is.matrix(ref_cells)) if (sum(is.na(ref_cells)) != 0) stop("'ref_cells' cannot hold NA values'")
    if (is.list(ref_cells)) if (sum(unlist(lapply(ref_cells, function (x) sum(is.na(x))))) != 0) stop("'ref_cells' cannot hold NA values'")

    cat('Preparing data....    \n')
    if (is.null(ref_cells)) {
        warning("Running without reference cells: classification only based on correlations \n", immediate. = TRUE)
    }

    ## If ref_cells is a list, make it one matrix and construct ref_types
    if (is.list(ref_cells)) {
        ref_types <- vector()
        for (i in seq_len(length(ref_cells))) {
            ref_types <- c(ref_types, rep(names(ref_cells)[i], ncol(ref_cells[[i]])))
        }
        ref_cells <- do.call(cbind, ref_cells)
        names(ref_types) <-colnames(ref_cells)
    }

    ## Same genes in reference as input
    if (!setequal(rownames(input), rownames(ref_cells))) {
        ref_cells <- equal.genes(input, ref_cells)
    }
    if (!setequal(rownames(input), rownames(ref_profiles))) {
       ref_profiles <- equal.genes(input, ref_profiles)
    }

    ## Make reference_profiles (one average profile per reference cell type):
    if (is.null(ref_profiles)) {
        ref_profiles <- MeanRef(ref_cells, method = "mean", ref_types = ref_types)
    }

    ## Make an environment to store the classification information and variables in
    Env <- new.env(parent = emptyenv())
    variables <- list(nodeheigths = c(), ##the heights of the different tree nodes, for plotting convenience
                    profilescores = list(), # a list of the profile scores per node
                    nodetypes = list(), # a list of the reference cell types per node
                    correlations = list(), # a list of the input correlations per node
                    confidencescores = list(), # a list of the confidence scores per node
                    genes = list(), # The gene names that are used in each node
                    counter = 0, # to keep track of the node depth/number
                    tree = NULL, # will be filled with the classification tree
                    p_thresh = p_thresh, # the p-value threshold for wilcox
                    fc_thresh = fc_thresh, # the fold-change threshold for wilcox (or fc, when fix_ngenes = FALSE)
                    pc_thresh = pc_thresh, # the % of expression threshold for wilcox
                    fix_ngenes = fix_ngenes, # boolean, use a fixed amount of genes?
                    subsample = subsample, # boolean, subsample ref_types, to min(table(ref_types))
                    only_pos = only_pos, # use only higher expressed genes
                    clust_dist = clust_dist, # clustering distance
                    cor_method = cor_method, # correlation method
                    gs_method = gs_method, # gene selection method (fold-change ('fc'), or wilcox)
                    print_steps = print_steps, # print the number of selected genes in each step?
                    clust_method = clust_method, # clustering method
                    n_genes = n_genes) # number of genes used, if fix_ngenes = TRUE
    Env <- list2env(variables, envir = Env)

    ## Calculate the profile and confidence scores in each node
    cat('Running analysis... \n')
    SplitNodes(Env = Env, # the environment
                input = input, # the input expression matrix
                all_ref_cells = ref_cells, ## matrix of all reference cells
                ref_types = ref_types, ## types of the all_ref_cells matrix
                ref_profiles =  ref_profiles) ## matrix of the average reference profiles per cell type

    ## For plotting purposes, find the x coordinates of the nodes of the classification tree
    heights <- Env$nodeheigths
    hc <- Env$tree
    coor <- cbind("y" = heights, "x" = rep(NA, length(heights)))
    coordinates <- dendextend::get_nodes_xy(as.dendrogram(hc))
    coordinates <- coordinates[coordinates[,2] != 0, ]
    for (i in seq_len(nrow(coor))) {
        x <- coordinates[ ,1][coordinates[ ,2] == coor[i,1]]
        coor[i,2] <- x
    }

    ## Define output
    parameters <- list(available.types = colnames(ref_profiles),
                        n_genes = n_genes,
                        cor_method = cor_method,
                        clust_method = clust_method,
                        clust_dist = clust_dist,
                        thresh = thresh,
                        subsample = subsample,
                        only_pos = only_pos,
                        pc_thresh = pc_thresh,
                        p_thresh = p_thresh,
                        fc_thresh = fc_thresh,
                        gs_method = gs_method,
                        fix_ngenes = fix_ngenes)

    output <- list(prof_scores = Env$profilescores,
                    conf_scores = Env$confidencescores,
                    nodetypes = Env$nodetypes,
                    correlations = Env$correlations,
                    tree = hc,
                    nodecoor = coor,
                    genes = Env$genes,
                    parameters = parameters)

    ## Finally, perform the actual classification
    output <- Classify(chetah = output,
                               thresh = thresh)

    if (plot.tree) PlotTree(chetah = output)
    return(output)
}     ### CHETAHclassifier
## -----------------------------------------------------------------------------

# Called by the CHETAHclassifier. Determines the branches
# of the current node of the classification tree,
# calls ScoreNode and reruns itself on the next-higher node(s).
# The fact that this function is recursive,
# makes it possible to do all calculations from top to bottom
SplitNodes <- function (input, Env, ref_profiles, all_ref_cells, ref_types) {

    ## Keep track of the number of the current node
    Env$counter <- Env$counter + 1

    ## (Re)construct the classification tree and cut at the highest node
    hc <- hclust(Env$clust_dist(t(ref_profiles)), method = Env$clust_method)
    branches <- cutree(hc, k = 2)

    ## Get the names of the celltypes in the two branches and save
    branch1 <- names(branches[branches == 1])
    branch2 <- names(branches[branches == 2])
    Env$nodetypes[[Env$counter]] <- branches
    Env$nodeheigths[Env$counter] <- max(hc$height)
    if (Env$counter == 1) Env$tree <- hc
    if (Env$print_steps) cat("Left node: ", branch1, "\nRight node: ", branch2, "\n")

    ## Calculate correlations, profile and confidence scores
    ScoreNode(branch1 = branch1,
                branch2 = branch2,
                input = input,
                Env = Env,
                ref_profiles = ref_profiles,
                all_ref_cells = all_ref_cells,
                ref_types = ref_types)

    ## If a branch contains 2 or more cell types (and thus other nodes),
    ## perform the steps in the lower nodes
    if (length(branch1) > 1) {
        SplitNodes(ref_profiles = ref_profiles[ ,branch1],input = input, Env = Env,
                   all_ref_cells = all_ref_cells, ref_types = ref_types)
    }
    if (length(branch2) > 1) {
        SplitNodes(ref_profiles = ref_profiles[ ,branch2], input = input, Env = Env,
                   all_ref_cells = all_ref_cells, ref_types = ref_types)
    }
}      ### SplitNodes
## -----------------------------------------------------------------------------

# Called by the CHETAHclassifier via SplitNode.
# Does the actual gene selections and correlations
# and calculates confidence scores.
ScoreNode <- function  (ref_profiles, branch1, branch2, input, Env, all_ref_cells, ref_types) {
    leaves <- c(branch1, branch2)
    Env$genes[[Env$counter]] <- list()

    ##################### First, gene selection and correlations are done for each reference cell type
    for(uniquetype in seq_len(length(leaves))) {

        ## save names of cell types in the same branch and cell types in the other branch
        current_type <- leaves[uniquetype]
        branch <- if(current_type %in% branch1) branch1 else branch2
        oth_branch <- if(current_type %in% branch1) branch2 else branch1

        ## Subsample types with many cells.
        if (!is.null(all_ref_cells)) {
            current_cells <- names(ref_types[ref_types %in% current_type])

            if (Env$subsample) {
                min_group_size <- min(table(ref_types[ref_types %in% oth_branch]))
                otherbranch_c <- c()
                for(i in seq_len(length(oth_branch))) {
                    cells <- names(ref_types[ref_types %in% oth_branch[i]]) ## ob == "other branch"
                    if (length(cells) > min_group_size ) {
                        cells <- cells[sample(length(cells), min_group_size)]
                    }
                    otherbranch_c <- c(otherbranch_c, cells)
                }
          } else {
              otherbranch_c <- names(ref_types[ref_types %in% oth_branch])
          }
        }

        ## select genes that are differentially expressed
        ## between current type and the other branch
        if (Env$gs_method == "wilcox" & !is.null(all_ref_cells)) {
            genes <- doWilcox(current_cells = current_cells,
                              otherbranch_c = otherbranch_c,
                              all_ref_cells = all_ref_cells,
                              Env = Env)
        } else {
            genes <- FindDiscrGenes(ref_profiles = ref_profiles,
                                    type = current_type,
                                    other_branch = oth_branch,
                                    Env = Env)
        }
        Env$genes[[Env$counter]][[current_type]] <- genes
        genes <- names(genes)

        ## Correlate to the reference cell type using the just selected genes
        if (Env$cor_method == "cosine") {
            correlations <- doCosine(genes = genes, input = input,
                                     all_ref_cells = all_ref_cells,
                                     current_cells = current_cells,
                                     otherbranch_c = otherbranch_c,
                                     ref_profiles = ref_profiles,
                                     type = current_type)
        } else {
            correlations <- doCorrelation(ref_profiles = ref_profiles,
                                          genes = genes, input = input,
                                          all_ref_cells = all_ref_cells,
                                          current_cells = current_cells,
                                          otherbranch_c = otherbranch_c,
                                          type = current_type,
                                          Env = Env)
        }
        ## correlation of input cells to the current reference profile (RP)
        cor_i <- correlations[["cor_i"]]
        ## correlation of reference cells of current type to RP
        cor_t <- correlations[["cor_t"]]
        ## correlation of reference cells in the other branch to RP
        cor_ob <- correlations[["cor_ob"]]

        ## produce the profile scores
        if(!is.null(all_ref_cells)) {
            good <- ecdf(cor_t)
            false <- ecdf(cor_ob)
            goodvalues <- good(cor_i) ; names(goodvalues) <- colnames(input)
            falsevalues <- (1 - false(cor_i)) ; names(falsevalues) <- colnames(input)
            profile_score <- goodvalues - falsevalues
        } else profile_score <- NA

        ## Add the profile scores and correlations of all types to one matrix
        if (uniquetype == 1) prof_scores <- as.matrix(profile_score) else {
            prof_scores <- cbind(prof_scores, profile_score)
        }
        if (uniquetype == 1) cors <- as.matrix(cor_i) else cors <- cbind(cors, cor_i)

    } ########### End of correlations per reference type

    if(Env$print_steps) cat("\n")

    colnames(cors) <- leaves
    if (!is.null(prof_scores)) colnames(prof_scores) <- leaves

    ## Calculate confidence scores
    if(!is.null(all_ref_cells)) {
      for(type_number in seq_len(length(leaves))) {

        ## Select type
        type <- leaves[type_number]

        ## Calculate the final confidence scores
        if(type %in% branch1) {
          prof_scores_ob <- apply(prof_scores[ ,branch2, drop = FALSE], 1, mean)
        } else {
          prof_scores_ob <- apply(prof_scores[ ,branch1, drop = FALSE], 1, mean)
        }
        confidence <- prof_scores[ ,type_number] - prof_scores_ob

        ## save confidence score
        if (type_number == 1) conf_scores <- as.matrix(confidence) else {
          conf_scores <- cbind(conf_scores, confidence)
        }
      }
      colnames(conf_scores) <- leaves
    } else conf_scores <- NA

    ## Save all the info of this step
    Env$confidencescores[[Env$counter]] <- conf_scores
    Env$correlations[[Env$counter]] <- cors
    Env$profilescores[[Env$counter]] <- prof_scores

}    ### ScoreNode
## -----------------------------------------------------------------------------
# Used by ScoreNode to do the gene selection, only based on fld-change between two groups of cells
FindDiscrGenes <- function (type, other_branch, ref_profiles, Env) {

    ## Calculate mean
    group1 <- apply(ref_profiles[ , type, drop = FALSE], 1, mean)
    group2 <- apply(ref_profiles[ , other_branch, drop = FALSE], 1, mean)
    group_diff <- group1-group2
    if (!Env$only_pos) {
        gd <- group_diff
        group_diff <- abs(group_diff)
    }

    ## select the genes
    if(!Env$fix_ngenes){
        pnn_genes <- group_diff
        pn_genes <- pnn_genes[pnn_genes > Env$fc_thresh]
        fc_thresh <- Env$fc_thresh
        while (length(pn_genes) == 0) {
            fc_thresh <- fc_thresh*0.5
            pn_genes <- pnn_genes[pnn_genes > fc_thresh]
        }
        pp_genes <- sum(pn_genes > 0)
        if (fc_thresh != Env$fc_thresh) cat("decreased fc_thresh to", fc_thresh, '   ')
        if(Env$print_steps) cat("|", length(pn_genes), "/", pp_genes, "|")
    } else {
        pn_genes <- sort(group_diff, decreasing = TRUE)[seq_len(Env$n_genes)]
        top <- sort(group_diff, decreasing = TRUE)[seq_len(Env$n_genes)]
        pn_genes <- rep(1, length(top))
        names(pn_genes) <- names(top)
        if (!Env$only_pos) pn_genes[gd[names(pn_genes)] < 0] <- 0
        # Print the amount of selected genes
        if(Env$print_steps) cat("|", length(pn_genes), "/", sum(pn_genes), "|")
    }

    return(pn_genes)
}        ### FindDiscrGenes

## -----------------------------------------------------------------------------
# Used by ScoreNode to do the gene selections using the Wilcox test
doWilcox <- function (current_cells, otherbranch_c, all_ref_cells, Env) {

    ## First calculate the fld change and percentage of cells expressing the gene, per gene.
    fld <- apply(all_ref_cells, 1, function (z) {
       mean(z[current_cells]) - mean(z[otherbranch_c])
    })
    pct1 <- apply(all_ref_cells, 1, function (z) {
        sum(z[current_cells] > 0) / length(current_cells)
    })
    pct2 <- apply(all_ref_cells, 1, function (z) {
        sum(z[otherbranch_c] > 0) / length(otherbranch_c)
    })

    ## Take only the genes in the all_ref_cells that exceed the fld and pct thresholds
    ## To prevent no genes exceeding the threshold, increase the threshold if necessary when needed
    fc_thresh <- Env$fc_thresh
    if(Env$only_pos) {
        testgenes  <- all_ref_cells[((pct1 > Env$pc_thresh | pct2 > Env$pc_thresh) & fld > fc_thresh), , drop = FALSE]
        while (nrow(testgenes) < 3) {
            fc_thresh <- 0.5 * fc_thresh
            testgenes <- all_ref_cells[((pct1 > Env$pc_thresh | pct2 > Env$pc_thresh) & (fld > fc_thresh)), , drop = FALSE]
        }
    } else {
        testgenes <- all_ref_cells[((pct1 > Env$pc_thresh | pct2 > Env$pc_thresh) & (fld > fc_thresh | fld < -fc_thresh )), , drop = FALSE]
        fc_thresh <- Env$fc_thresh
        while (nrow(testgenes) < 3) {
            fc_thresh <- 0.5 * fc_thresh
            testgenes <- all_ref_cells[((pct1 > Env$pc_thresh | pct2 > Env$pc_thresh) & (fld > fc_thresh | fld < -fc_thresh )), , drop = FALSE]
        }
    }
    if (fc_thresh != Env$fc_thresh) cat("decreased fc_thresh to", log2(fc_thresh), '   ')

    ## Wilcox test + p-adjustment & select the genes
    pval <- apply(testgenes, 1, function (z) {
        wtest <- wilcox.test(x = z[current_cells], y = z[otherbranch_c])
        wtest$p.value
    })
    pval[is.na(pval)] <- 1
    pval <- p.adjust(pval, method = "bonferroni", n=length(pval))

    genes <- names(pval)[pval < Env$p_thresh]

    ## To prevent no genes exceeding the threshold, increase the threshold if necessary when needed
    p_thresh <- Env$p_thresh
    while (length(genes) == 0) {
        p_thresh <- p_thresh * 10
        genes <- names(pval)[pval < p_thresh]
    }
    if (p_thresh != Env$p_thresh) cat("p-val increased to", p_thresh, '   ')

    ## which genes are higher expressed and which lower
    pn_genes <- rep(0, length(genes))
    names(pn_genes) <- genes
    pn_genes[fld[genes] > 0] <- 1

    ## Print the amount of selected genes
    if (Env$print_steps) cat("|", length(pn_genes), "/", sum(pn_genes), "|")

    return(pn_genes)
}    ## doWilcox
## -----------------------------------------------------------------------------
# Used by ScoreNode to do the correlations using cosine
doCosine <- function(genes, input, all_ref_cells, current_cells,
                     otherbranch_c, ref_profiles, type) {

    ## matrices with selected genes + cells
    cells_i <- input[genes, , drop = FALSE]
    cells_t <- ref_profiles[genes, type, drop = FALSE]

    ## Do cosine
    t_vs_t <- drop(crossprod(cells_t,cells_t))
    cor_i <-  drop(crossprod(cells_t,cells_i)/sqrt(diag(crossprod(cells_i,cells_i)) * t_vs_t))
    cor_i[is.na(cor_i)] <- 0

    names(cor_i) <- colnames(input)

    if(!is.null(all_ref_cells)) {
        ## same steps for type and other branch
        cells_b <- all_ref_cells[genes , current_cells, drop = FALSE]
        cells_ob <- all_ref_cells[genes , otherbranch_c, drop = FALSE]
        cor_t <-  drop(crossprod(cells_t,cells_b)/
                         sqrt(diag(crossprod(cells_b,cells_b)) * t_vs_t))
        cor_ob <- drop(crossprod(cells_t,cells_ob)/
                         sqrt(diag(crossprod(cells_ob,cells_ob)) * t_vs_t))
        cor_t[is.na(cor_t)] <- 0
        cor_ob[is.na(cor_ob)] <- 0
        names(cor_t) <- current_cells
        names(cor_ob) <- otherbranch_c
    } else {
        cor_t <- NULL
        cor_ob <- NULL
    }

    data <- list(cor_i, cor_t, cor_ob)
    names(data) <- c("cor_i", "cor_t", "cor_ob")
    return(data)
}     ##doCosine
## -----------------------------------------------------------------------------
# Used by ScoreNode to do the correlations (if Env$cor_method != cosine)
doCorrelation <- function(ref_profiles, genes, input, all_ref_cells,
                          current_cells, otherbranch_c, type, Env) {
    ## as.matrix: in cases were a sparse Matrix is used
    ## suppressWarnings: suppress warnings when the sd is 0 and cor returns NA
    cor_i <-   suppressWarnings(cor(as.matrix(input[genes, , drop = FALSE]),
                                    as.matrix(ref_profiles[genes, type, drop = FALSE]),
                                    method = Env$cor_method))
    cor_i[is.na(cor_i)] <- 0
    names(cor_i) <- colnames(input)

    if(!is.null(all_ref_cells)) {
        cor_t <- suppressWarnings(cor(as.matrix(all_ref_cells[genes , current_cells, drop = FALSE]),
                                      as.matrix(ref_profiles[genes, type, drop = FALSE]),
                                      method = Env$cor_method))
        cor_ob <- suppressWarnings(cor(as.matrix(all_ref_cells[genes , otherbranch_c, drop = FALSE]),
                                       as.matrix(ref_profiles[genes, type, drop = FALSE]),
                                       method = Env$cor_method))
        cor_t[is.na(cor_t)] <- 0
        cor_ob[is.na(cor_ob)] <- 0
        names(cor_t) <- current_cells
        names(cor_ob) <- otherbranch_c
    } else {
        cor_t <- NULL
        cor_ob <- NULL
    }

    data <- list(cor_i, cor_t, cor_ob)
    names(data) <- c("cor_i", "cor_t", "cor_ob")
    return(data)
}     ## doCorrelation
## -----------------------------------------------------------------------------
#' (Re)classify after running \code{\link{CHETAHclassifier}} using a confidence threshold \cr
#' NOTE: In case of bulk reference profiles: only the correlations will be used,
#' as the data does not allow for profile or confidence scores to be calculated.
#'
#' @param chetah the CHETAH object
#' @param thresh a confidence threshold between -0 and 2. \cr
#' Selecting 0 will classify all cells, whereas 2 will result i
#' n (almost) no cells to be classified. \cr
#' \emph{recommended}: between 0.1 (fairly confident) and 1 (very confident)
#' @return_clas Instead of returning the chetah object, only return the classification vector
#' @return
#' a charachter vector of the cell types with the names of the cells
#' @export
#'
#' @examples
#' ## Classify all cells
#' chetah <- Classify(chetah, 0)
#'
#' ## Classify only cells with a very high confidence
#' chetah <- Classify(chetah, 1)
Classify <- function(chetah, thresh, return_clas = FALSE) {

    ## Get parameters
    if(!is.na(chetah$conf_scores[[1]][1])) {
        conf <- chetah$conf_scores
        prof <- chetah$prof_scores
    } else {
        conf <- NULL
        prof <- chetah$correlations
    }

    nodetypes <- chetah$nodetypes
    classification <- nodeDown(conf = conf, prof = prof, node = 1,
                               thresh = thresh, nodetypes = nodetypes)
    names(classification) <- rownames(prof[[1]])
    chetah$classification <- classification
    if (return_clas) return(classification) else return(chetah)
}

## ---------------------------------
# To make a profile_matrix from a list of reference matrices
MeanRef <- function (input,
                     method = "mean",
                     scale = 10000,
                     log = FALSE,
                     ref_types) {

    if (!(method == "mean" | method == "median" | method == "var")) {
        stop ("method must be either mean or median")
    }
    if (is.list(input)) {
        ## log normalize
        if(log) {
            ref_cells <- lapply(input, function(x) apply(x, 2, function(column) log2( (column/sum(column) * scale) +1 )))
        } else ref_cells <- input

        ## Average each cell type
        ref_cells <- lapply( input, function (x) apply(x, 1, method) )
        ref_means <- do.call(cbind, ref_cells)

        ## Take the cell types as colnames and the genes as rownames
        colnames(ref_means) <- names(input)
        rownames(ref_means) <- rownames(input[[1]])
      } else {
            ref_means <- matrix(nrow = nrow(input))
            for(i in unique(ref_types)) {
                new <- apply(input[ ,ref_types == i, drop = FALSE], 1, method)
                ref_means <- cbind(ref_means, new)
            }
            ref_means <- ref_means[ ,-1]
            colnames(ref_means) <- unique(ref_types)
            rownames(ref_means) <- rownames(input)
    }
    ## return
    return(ref_means)
}              ## MeanRef
## ---------------------------------------------------------------------------------------------------------------------------

#' Plots the chetah classification tree with nodes numbered
#'
#' @param chetah a chetah object
#' @param col a vector of colors, with the names of the reference cell types
#' @param col_nodes a vector of colors, ordered for node 1 till the last node
#' @param return instead of printing, return the ggplot object
#' @param no_bgc remove the background color from the node numbers
#' @param plot_limits define the % of the plot_heigth that should be used as a margin below and above the plot
#' Decreasing the former further is usefull when the labels are cut of the plot (default = c(-0,25, 01)).
#' @param labelsize the size of the intermediate and leaf node labels (default = 6)
#' @return
#' A ggplot object of the classification tree
#' @export
#' @import ggplot2
#' @importFrom dendextend color_branches
#' @importFrom dendextend color_labels
#' @importFrom stats as.dendrogram
#' @examples PlotTree(chetah = chetah)
PlotTree <- function(chetah, col = NULL,
                     col_nodes = NULL,
                     return = FALSE,
                     no_bgc = FALSE,
                     plot_limits = c(-0.4, 0.1),
                     labelsize = 6) {
    ## Extract parameters + make dendrogram
    lbls <- chetah$tree$labels
    ord <- chetah$tree$order
    chetah$tree$labels <- paste0(chetah$tree$labels, "  ")
    hc <- as.dendrogram(chetah$tree)
    hc <- dendextend::set(hc, 'labels_cex', labelsize/6)

    ## Give colors
    if(!is.null(col)) {
        col1 <- col[lbls[ord]]
        hc <- dendextend::color_branches(hc, clusters = seq_len(length(lbls)), col = col1)
        hc <- dendextend::color_labels(hc, col = col1)
    }
    if(!is.null(col_nodes)) {
        col_n <- col_nodes[grepl("Node", names(col_nodes))]
        names(col_n) <- gsub("Node", "", names(col_n))
        col_n <- col_n[order(as.numeric(names(col_n)))][seq_len(nrow(chetah$nodecoor))]
        col_r <- col_n[seq_len(nrow(chetah$nodecoor))]
        alpha1 <- 0.2
    } else {
        col_n <- rep("black", nrow(chetah$nodecoor))
        col_r <- "white"
        alpha1 <- 0.8
    }
    if(no_bgc) {
        col_r <- "white"
        alpha1 <- 0.8
    }

    ## Extract additional parameters
    maxh <- max(dendextend::get_branches_heights(hc))
    length <- length(labels(hc))

    ## Plot
    ggdend <- ggplot(hc) +
      annotate(geom = "rect",
               xmin = chetah$nodecoor[, "x"] - 0.02*length,
               xmax = chetah$nodecoor[, "x"] + 0.02*length,
               ymin = chetah$nodecoor[, "y"] - 0.05*maxh,
               ymax = chetah$nodecoor[, "y"] + 0.05*maxh,
               fill = col_r, alpha = alpha1) +
      annotate(geom = "text",
               x = chetah$nodecoor[, "x"],
               y = chetah$nodecoor[, "y"],
               label = seq_len(nrow(chetah$nodecoor)),
               size = labelsize, fontface = "bold",
               col = col_n) +
      ylim(plot_limits[1]*maxh, maxh+plot_limits[2]*maxh) +
      ggtitle("Classification Tree")

    ## Delete layer that gives warnings
    keep <- !unlist(lapply(ggdend$layers, function(x) class(x$geom)[1] == "GeomPoint"))
    ggdend$layers <- ggdend$layers[keep]
    if (return) return(ggdend) else suppressWarnings(print(ggdend))
}   ## PlotTree
## -----------------------------------------------------------
#' Plots a variable on a t-SNE
#'
#' @param toplot the variable that should be plotted. Either a character vector
#' (which will become a factor) or a factor, or a (continuous) numeric.
#' If toplot is not named with the rownames of \code{coor}, it is assumed
#' that the order of the two is the same.
#' @param coor the tSNE, or other 2D vizualization coordinates.
#' A matrix or dataframe with cells in the rows and x (e.g. tSNE_1) and
#' y (e.g. tSNE_2) coordinates in the 1st and 2nd column respectively
#' @param col a vector of colors. If \code{toplot} is a numeric,
#' this will become a continuous scale. \cr
#' \emph{If \code{toplot} is a charachter vector, the colors should
#' be named with the unique values (/levels) of toplot}
#' @param return instead of printing, return the ggplot object
#' @param limits the limits of the continuous variable to plot.
#' When not provided the minimal and maximal value will be used
#' @param pt.size the point-size
#' @param shiny Needed for the shiny application: should always be NULL
#' @param y_limits the y-axis limits
#' @param x_limits the x-axis limits, if NULL
#' @param legend_label the label of the legend
#' @return
#' A ggplot object
#' @export
#' @import ggplot2
#' @importFrom gplots colorpanel
#' @examples
#' tsne <- tsne_mel
#' CD8 <- as.matrix(data_mel['CD8A', ])
#' PlotTSNE(toplot = CD8, coor = tsne)
PlotTSNE <- function (toplot, coor, col = NULL, return = FALSE,
                      limits = NULL, pt.size = 1, shiny = NULL,
                      y_limits = NULL, x_limits = NULL, legend_label = '') {
    ## Merge coor and toplot
    if (is.null(names(toplot)) & is.null(dim(toplot))) names(toplot) <- rownames(coor)
    toplot <- data.frame(toplot, stringsAsFactors = TRUE)
    if (!is.null(shiny) & ncol(toplot) == 2) colnames(toplot)[2] <- 'key'
    toplot$rn <- rownames(toplot)
    colnames(toplot)[1] <- "variable"
    coor <- data.frame(coor)
    coor$rn <- rownames(coor)
    colnames(coor) <- c("tSNE_1", "tSNE_2", "rn")
    data <- merge(toplot, coor, by = "rn")

    ## Initial plot
    if (!is.null(shiny)) {
        text <- paste0('cell_label: ', data$rn, '<br>', data$shiny, data$variable)
        if ('key' %in% colnames(data)) text <- paste0(text, '<br>','cell type: ', data$key)
        data$text <- text
        plot <- ggplot(data, aes_string(x='tSNE_1', y='tSNE_2', text = 'text'))
    } else {
        plot <- ggplot(data, aes_string(x='tSNE_1', y='tSNE_2'))
    }
    plot <- plot +
        theme_classic() +
        geom_point(aes_string(color = 'variable'), size = pt.size) +
        labs(color = legend_label) +
        theme(legend.text = element_text(size=10),
              legend.title = element_text(size=10))

    class <- class(toplot[,1])

    # For a categorical variable:
    if(class == "factor") {
      plot <- plot +
        guides(colour = guide_legend(override.aes = list(size = 2.5*pt.size),
                  ncol = 1, title = legend_label))
      if (!is.null(col)) plot <- plot + scale_color_manual(values = col)
    }

    # For a numeric variable:
    if(class == "numeric" | class == 'integer') {
      if (is.null(limits)) {
        limits <- c(min(toplot[,1])-(0.01*min(toplot[,1])),
                      max(toplot[,1])+(0.01*max(toplot[,1])))
      }
      if (is.null(col)) {
        col <- c(gplots::colorpanel(n = 50, low = '#2f2bad', mid = '#93cfe2', high = '#fffaaa'),
                  gplots::colorpanel(n = 50, low = '#fffaaa', mid = "#ffad66", high = '#d60000'))
      }
      plot <- plot +
          scale_color_gradientn(colors = col, limits = limits)
    }

    ## Define plotting limits
    if (!is.null(x_limits)) {
        plot <- plot + lims(x = c(x_limits[1]-1, x_limits[2]+1)) }
    if (!is.null(y_limits)) {
        plot <- plot + lims(y = c(y_limits[1]-1, y_limits[2]+1)) }
    if (return) return(plot) else print(plot)
}
## ----------------------------------------------------------------------------
#' Plot boxplots grouped by a class variable
#'
#' @param toplot a numeric variable. When provided as a dataframe,
#'  the colname will be used a a title
#' @param class a factor of the length of \code{toplot}
#' to indicate the groups
#' @param col a character vector of colors, with the
#' names of the factor levels
#' @param grad_col if defined, the jittered points of
#' the individual values of each cell behind
#' the boxplot will be colored by these gradient colors
#' based on their value
#' @param limits the limits of \code{toplot}
#' @param return instead of printing, return the ggplot object
#'
#' @return
#' A ggplot object
#' @import ggplot2
#' @importFrom grDevices rgb
#' @export
#' @examples
#' CD8 <- as.matrix(data_mel['CD8A', ])
#' types <- chetah$classification
#' PlotBox(toplot = CD8, class = types)
PlotBox <- function(toplot, class, col = NULL, grad_col = NULL,
                    limits = NULL, return = FALSE) {

    toplot <- data.frame(toplot)
    class <- data.frame(class, stringsAsFactors = TRUE)
    col <- col[levels(class[,1])]

    if (is.null(limits)) {
        limits <- c(min(toplot[,1])-0.01, max(toplot[,1])+0.01)
    }
    data <- cbind.data.frame(toplot, class)
    colnames(data) <- c("score", "class")

    ## Plot basics + background
    plot <- ggplot(data, aes_string(x = 'class', y = 'score')) +
        geom_hline(yintercept = 0, color = "gray80") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10)) +
        labs(color = "cell types") +
        xlab("cell types") +
        ylab("score") +
        ylim(low = limits[1], high = limits[2]) +
        coord_flip()

    ## background jitter (value per cell, scattered on the x-axis)
    if (!is.null(col)) plot <- plot + scale_color_manual(values = col)
    if (is.null(grad_col)) {
        plot <- plot + geom_jitter(aes_string(color = 'class'), size = 0.1)
    } else {
        plot <- plot + geom_jitter(aes_string(fill = 'score'),
                                    size = 0.1, shape = 21, color = rgb(0,0,0,0)) +
        scale_fill_gradientn(colors = grad_col, limits = c(-2,2))
    }
    ## Add the boxplot
    plot <- plot  + geom_boxplot(aes_string(color = 'class'),
                                 fill = rgb(1,1,1, 0.6), outlier.colour="NA")

    if (return) return(plot) else print(plot)
}
## -------------------------------------------------------------------------------------
#' Plot the CHETAH classification on 2D visulization like t-SNE
#' + the corresponding classification tree,
#' colored with the same colors
#'
#' @param chetah a chetah object
#' @param coor the tSNE, or other 2D vizualization coordinates.
#' A matrix or dataframe with cells in the rows and x (e.g. tSNE_1)
#'  and y (e.g. tSNE_2) coordinates in the 1st and 2nd column respectively
#' @param interm color the intermediate instead of the final types
#' @param tree plot the tree, along with the classification
#' @param pt.size the point-size of the classication plot
#' @param return_col whether the colors that are used for
#' the classification plot should be returned
#' @param col custom colors for the cell types. \emph{the colors
#' should be named with the corresponding cell types}
#' @param return return the plot instead of printing it
#' @return a ggplot object
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom grDevices rainbow
#' @examples
#' tsne <- tsne_mel
#' ## Standard plot (final types colored)
#' PlotCHETAH(chetah = chetah, coor = tsne)
#'
#' ## Intermediate types colored
#' PlotCHETAH(chetah = chetah, coor = tsne, interm = TRUE)
#'
#' ## Plot only the t-SNE plot
#' PlotCHETAH(chetah = chetah, coor = tsne, tree = FALSE)
PlotCHETAH <- function(chetah, coor, interm = FALSE, return = FALSE,
                       tree = TRUE, pt.size = 1, return_col = FALSE, col = NULL) {
    ## Determine colors: if not custom, use the predefined ones. If too many cell types: rainbow.
    if(is.null(col)) {
        colors <- c ('blue', 'gold', 'cyan3', 'navy',
                   'forestgreen', 'orange', 'darkolivegreen3',
                   'brown', 'green', 'purple','deepskyblue', 'cyan',
                   'orangered3', 'coral', 'yellow3', "black",
                   'yellow1', 'darkorchid1', 'darksalmon', 'darkseagreen1',
                   'darkslategrey', 'deeppink4', 'green2', 'lemonchiffon1',
                   'lightcyan', 'midnightblue', 'maroon1', 'orange3', 'palegreen',
                   'palevioletred1', 'peru', 'seagreen1', 'red3', 'snow2',
                   'steelblue1', 'turquoise')
        leaf_nodes <- c(names(chetah$nodetypes[[1]]), unique(chetah$classification)[!(unique(chetah$classification) %in% chetah$nodetypes[[1]])])
        int_nodes <- paste0("Node", seq_len(length(chetah$nodetypes)))
        M <- max(length(int_nodes), length(leaf_nodes))
        if (length(int_nodes) < 24)   {
            gray <- paste0('gray', rev(seq(2, 92, 4)))
        } else if (length(int_nodes) < 47) {
            gray <- paste0('gray', rev(seq(2, 92, 2)))
        } else {
            gray <- rep('gray70', M)
        }
        if (M > 36) colors <- grDevices::rainbow(M)

        # color either the final, or intermediate types
        if(!interm) {
            names(colors) <- leaf_nodes
            names(gray) <- int_nodes
        } else {
            names(gray) <- leaf_nodes
            names(colors) <- int_nodes
        }
        col <- c(gray, colors)
        col <- col[!is.na(names(col))]
    } else col <- col

    ## Plot
    plot1 <- PlotTSNE(toplot = chetah$classification, coor = coor,
                    col = col, return = TRUE, pt.size = pt.size)
    if(!interm) plot2 <- PlotTree(chetah, col, return = TRUE)
    if(interm) plot2 <- PlotTree(chetah, col_nodes = col, no_bgc = TRUE, return = TRUE)
    if (tree) {
        plots <- cowplot::plot_grid(plot1, plot2, ncol = 2)
    } else plots <- plot1
    if (return) return(plots) else print(plots)
    if (return_col) return(col)
}
## ----------------------------------------------------------------------------------------
#' Correlate all reference profiles to each other
#' using differentially expressed genes.
#'
#' @param ref_cells the reference, similar to
#' \code{\link{CHETAHclassifier}}'s ref_cells
#' @param ref_types must be provided if \code{ref_cells} is one matrix
#' @param ref_profiles similar to
#' \code{\link{CHETAHclassifier}}'s ref_profiles
#' @param return return the matrix that was used to produce the plot
#' @param n_genes as in \code{\link{CHETAHclassifier}}
#' @param fix_ngenes as in \code{\link{CHETAHclassifier}}
#' @param print_steps as in \code{\link{CHETAHclassifier}}
#' @param only_pos as in \code{\link{CHETAHclassifier}}
#' @return
#' A square plot. The values show how much two reference profiles
#' correlate, when using the genes with the highest fold-change.
#'
#' @export
#'
#' @importFrom gplots colorpanel
#' @importFrom corrplot corrplot
#' @importFrom stats cor
#' @examples
#' reference <- reference_hn
#' CorrelateReference(ref_cells = reference)
CorrelateReference <- function(ref_cells = NULL, ref_types = NULL,
                                ref_profiles = NULL, return = FALSE, n_genes = 200,
                                fix_ngenes = TRUE, print_steps = FALSE, only_pos = FALSE) {
    cat('Running... this may take a couple of minutes \n')
    stopifnot(is.list(ref_cells) | (!is.null(ref_types) &
                !is.null(ref_cells)) | is.null(ref_profiles))

    ## Make ref_profiles
    if (is.null(ref_profiles)) ref_profiles <- MeanRef(ref_cells, ref_types = ref_types, method = 'mean')

    ## Make empthy correlation matrix
    cors <- matrix(NA, nrow = ncol(ref_profiles), ncol = ncol(ref_profiles))
    rownames(cors) <- colnames(ref_profiles)
    colnames(cors) <- colnames(ref_profiles)
    Env <- list('n_genes' = n_genes, 'only_pos' = only_pos,
                'print_steps' = print_steps, 'fix_ngenes' = fix_ngenes)

    ## Do the pair-wise correlations
    for(i in seq_len(ncol(ref_profiles))) {
        for(j in seq_len(ncol(ref_profiles))) {
            if (i == j) cors[i,j] <- 1 else {
                if (!is.na(cors[j,i])) cors[i,j] <- cors[j,i] else {
                genes <- names(FindDiscrGenes(type = colnames(ref_profiles)[i],
                                                other_branch = colnames(ref_profiles)[j],
                                                ref_profiles = ref_profiles,
                                                Env = Env))

                corr <- cor(ref_profiles[genes,i, drop = FALSE],
                            ref_profiles[genes, j, drop = FALSE],
                            method = "spearman")
                cors[i,j] <- corr
                }
            }
        }
    }
    ## plot
    cols <- c(gplots::colorpanel(n = 10, low = '#2f2bad', mid = '#93cfe2', high = '#fffdd3'),
              gplots::colorpanel(n = 10, low = '#fffdd3', mid = "#ffad66", high = '#d60000'))
    corrplot::corrplot(cors, method = "square", order = "hclust",
                       col = cols, mar = c(0, 0, 2, 0), type = 'upper')

    if (return) return(cors)
} ## CorrelateReference
## ----------------------------------------------------------------------------------------
#' Use a reference dataset to classify itself.
#' A good reference should have almost no mixture
#' between reference cells.
#'
#' @param ref_cells the reference, similar to
#' \code{\link{CHETAHclassifier}}'s ref_cells
#' @param ref_types must be provide if
#' \code{ref_cells} is a matrix of dataframe
#' @param return return the matrix that was used to produce the plot
#' @param ... Other variables to pass to
#' \code{\link{CHETAHclassifier}}
#' @return
#' A square plot. The rows are the original cell types,
#' the columns the classifion labels.
#' The colors and sizes of the squares indicate which
#' part of the cells of the rowname type are
#' classified to the type of the column name.
#' On the left of the plot, the percentage of cells
#' that is classified to an intermediate type
#' is plotted.
#' A good reference would classify nearly 100% of cells of type A to type A.
#' @export
#'
#' @importFrom gplots colorpanel
#' @importFrom corrplot corrplot
#' @importFrom graphics text
#' @examples
#' reference <- reference_hn
#' ClassifyReference(ref_cells = reference)
ClassifyReference <- function(ref_cells, ref_types = NULL, return = FALSE, ...) {
    stopifnot(is.list(ref_cells) | !is.null(ref_types))

    ## Prepare the reference
    if(is.list(ref_cells)) {
        ref_types <- c(rep( names(ref_cells)[1], ncol(ref_cells[[1]])),
                       rep( names(ref_cells)[2], ncol(ref_cells[[2]])))
        for (i in 3:length(ref_cells)) {
            ref_types <- c(ref_types, rep(names(ref_cells)[i], ncol(ref_cells[[i]])))
        }
        ref_cells <- do.call(cbind, ref_cells)
        names(ref_types) <-colnames(ref_cells)
    }

    ## Classify
    chetah <- CHETAHclassifier(input = ref_cells, ref_cells = ref_cells,
                                ref_types = ref_types, ...)

    ## Make the correlation matrix
    nms <- unique(ref_types)
    lngt <- length(nms)
    type <- chetah$classification
    splt <- unique(chetah$classification)
    splt <- splt[grepl("Node", splt)]
    cors <- matrix(NA, nrow = lngt, ncol = lngt+1)
    rownames(cors) <- nms
    colnames(cors) <- c(nms, 'Nodes')

    ## Extract the percentages
    for(i in seq_len(lngt)) {
        for(j in seq_len(lngt)) {
            actual <- names(ref_types)[ref_types == nms[i]]
            cors[i,j] <- sum(type[actual] == nms[j])/length(actual)
        }
    }

    ## Extract the percentage that ended up in a node
    for(i in seq_len(lngt)) {
        actual <- names(ref_types)[ref_types == nms[i]]
        cors[i,lngt+1] <- sum(grepl("Node", type[actual]))/length(actual)
    }

    ## Plot
    cols <- c(rep('black', 20),
                    gplots::colorpanel(n = 10, low = '#2f2bad', mid = '#93cfe2', high = '#fffdd3'),
                    gplots::colorpanel(n = 10, low = '#fff6b5', mid = "#ffad66", high = '#d60000'))
    a <- corrplot::corrplot(cors[ ,seq_len(lngt)], method = "square", order = "hclust",
                            col = cols, cl.lim = c(0,1), mar = c(0, 0, 2, 0))
    text(-3.5, seq_len(nrow(cors)), labels = round(cors[ ,'Nodes'], 2)[rev(rownames(a))])
    text(-3.5, nrow(cors)+1, '(%) in nodes')

    if (return) return(cors)
} ## ClassifyReference
## ----------------------------------------------------------------------------------
nodeDown <- function(conf, prof, node, thresh, nodetypes, prev_clas = NULL) {

    ## Do the classification of this node
    sub_clas <- apply(prof[[node]], 1, which.max)
    sub_clas  <- colnames(prof[[node]])[sub_clas]
    if (!is.null(conf)) {
        sub_score <- c()
        for (row in seq_len(nrow(conf[[node]]))) {
          sub_score <- c(sub_score, conf[[node]][row, sub_clas[row]])
        }
        sub_clas[sub_score < thresh] <- paste0("Node", node)
    }
    ## Only apply on cells reaching this node
    if(!is.null(prev_clas)) {
        reclas <- prev_clas %in% names(nodetypes[[node]])
        sub_clas[!reclas] <- prev_clas[!reclas]
    }
    ## continue to next node, or return
    if(node != length(prof)){
        nodeDown(conf = conf,
                 node = node + 1,
                 thresh = thresh,
                 nodetypes = nodetypes,
                 prof = prof,
                 prev_clas = sub_clas)
    } else {
        return(sub_clas)
    }
}
## select the common genes of the input and the reference
equal.genes <- function(inp, ref) {
    common <- intersect(rownames(inp), rownames(ref))
    ref <- ref[common, ]
}
### ---------------------------------------------------------------------------
## Select all the final AND intermediate types below a node
SelectNodeTypes <- function (chetah, whichnode) {
  nodes <- unlist(lapply(chetah$nodetypes[(whichnode):length(chetah$nodetypes)], function (x) {
    all(names(x) %in% names(chetah$nodetypes[[whichnode]]))
  }))
  nodes <- paste0('Node', ((whichnode):length(chetah$nodetypes))[nodes])
  alltypes <- c(nodes, names(chetah$nodetypes[[whichnode]]))
  return(alltypes)
}
### ---------------------------------------------------------------------------
#' In the CHETAH classification, replace the name of a Node
#' and all the names of the final and intermediate types under that Node.
#'
#' @param chetah a chetah object
#' @param whichnode the number of the Node
#' @param replacement a character vector that replaces the names under the selected Node
#' @param nodes_exclude \emph{optional} numbers of the Nodes under the selected Node, that should \strong{NOT} be replaced
#' @param node_only only rename the Node itself, without affecting the types under that Node
#' @param return_clas Instead of returning the chetah object, only return the classification vector
#'
#' @return
#' The chetah object with the new classification or if `return_clas = TRUE` the classification vector.
#' @export
#'
#' @examples
#' In the example data replace all T-cell subtypes by "T cell"
#' chetah <- RenameBelowNode(chetah = chetah, whichnode = 7, replacement = "T cell")
RenameBelowNode <- function(chetah, whichnode, replacement, nodes_exclude = NULL, node_only = FALSE, return_clas = FALSE) {
    classification <- chetah$classification
    if (node_only) {
        classification[classification == paste0("Node", whichnode)] <- replacement
    } else {
        replace <- SelectNodeTypes(chetah = chetah, whichnode = whichnode)
        if (!is.null(nodes_exclude)) {
            exclude <- vector()
            for (number in 1:length(nodes_exclude)) {
                exclude <- c(exclude, SelectNodeTypes(chetah = chetah, whichnode = nodes_exclude[number]))
            }
        replace <- replace[!(replace %in% exclude)]
        }
        classification[classification %in% replace] <- replacement
        chetah$classification <- classification
        if (return_clas) return(classification) else return(chetah)
    }
}
### -----------------
ch_env <- new.env(parent = emptyenv())
### ---------------------------------------------------------------------------
## Package environment to transfer the data for the shiny package
#' Launch a web page to interactively go trough the classification
#'
#' @param chetah a chetah object
#' @param coor the coordinates of 2D visualization,
#' e.g. t-SNE. 2 columns for the 2 coordinates types
#' (e.g. tSNE_1, tSNE_2). The column names have
#' to overlap with the names of the input count data of chetah.
#' @param counts the counts of the data that was
#' the input to the chetah classification
#' @importFrom plotly renderPlotly plotlyOutput ggplotly event_data
#' @import shiny
#' @importFrom pheatmap pheatmap
#' @importFrom reshape2 melt
#' @return
#' Opens a web page in your default browser
#' @export
#'
#' @examples
#' \dontrun{
#' counts <- data_mel
#' tsne <- tsne_mel
#' CHETAHshiny(chetah = chetah, coor = tsne, counts = counts)
#' }
CHETAHshiny <- function (chetah, coor, counts) {
    ## Search for the shinyApp in the package directory
    appDir <- system.file("shinyApp", "App", package = "CHETAH")
    if (appDir == "") {
      stop("Could not find example directory. Try re-installing `chetah`.", call. = FALSE)
    }
    ## Put data in package_environment
    ch_env$chetah <- chetah
    ch_env$coor <- coor
    ch_env$counts <- counts
    rm(chetah, coor, counts)
    shiny::runApp(appDir, display.mode = "normal")
}

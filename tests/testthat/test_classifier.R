context("Classifier")

chetah <- CHETAHclassifier(input = input_mel, ref_cells = headneck_ref)

test_that("Output classes", {
    expect_is(chetah@int_metadata$CHETAH$nodetypes, "list")
    expect_is(chetah@int_metadata$CHETAH$tree, "hclust")
    expect_is(chetah@int_metadata$CHETAH$nodecoor, "matrix")
    expect_is(chetah@int_metadata$CHETAH$nodecoor, "matrix")
    expect_is(chetah@int_metadata$CHETAH$genes, "list")
    expect_is(chetah@int_metadata$CHETAH$parameters, "list")
    expect_is(chetah@int_metadata$CHETAH$input_c, "character")
    expect_is(chetah@int_colData$CHETAH$prof_scores, "DataFrame")
    expect_is(chetah@int_colData$CHETAH$conf_scores, "DataFrame")
    expect_is(chetah@int_colData$CHETAH$correlations, "DataFrame")
    expect_is(chetah$celltype_CHETAH, "character")
})

test_that("Output consistency", {
    expect_equal(colnames(chetah@int_colData$CHETAH$correlations), colnames(chetah@int_colData$CHETAH$prof_scores))
    expect_equal(colnames(chetah@int_colData$CHETAH$prof_scores), colnames(chetah@int_colData$CHETAH$prof_scores))
    expect_equal(nrow(chetah@int_colData$CHETAH$prof_scores), ncol(chetah))
    expect_equal(names(chetah$celltype_CHETAH), colnames(chetah))
    expect_true(all(chetah$celltype_CHETAH %in% c(names(chetah@int_metadata$CHETAH$nodetypes[[1]]), 
                                      colnames(chetah@int_colData$CHETAH$prof_scores), 
                                      "Unassigned")))
})





    
context("Plots")
test_that("return colors", { 
    expect_is(PlotCHETAH(input_mel, return_col = TRUE), "character")
})
test_that("plotchetah output", {
    vdiffr::expect_doppelganger("plot: standard chetah", PlotCHETAH(input_mel))
    vdiffr::expect_doppelganger("plot: intermediate types", PlotCHETAH(input = input_mel, interm = TRUE))
    vdiffr::expect_doppelganger("plot: tree only", PlotCHETAH(input = input_mel, tree = FALSE))
    vdiffr::expect_doppelganger("plot: rename below", {
        chetah_rb <- RenameBelowNode(input = input_mel, whichnode = 6, replacement = "TCELL")
        PlotCHETAH(input = chetah_rb, tree = FALSE)
    })
    vdiffr::expect_doppelganger("plot: reclassify", {
        chetah_reclass <- Classify(input = input_mel, thresh = 0.5)
        PlotCHETAH(input = chetah_reclass, tree = FALSE)
    })
})
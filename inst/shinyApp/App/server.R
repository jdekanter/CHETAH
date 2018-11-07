server <- function(input, output, session) {
  ## ----------------------------------- Select colors ------------------
  counts <- CHETAH:::ch_env$counts
  chetah <- CHETAH:::ch_env$chetah
  coor   <- CHETAH:::ch_env$coor
  rm(counts, chetah, coor, envir = ch_env)
  gc()
  ## Confidence t-SNE color gradient
  grad_col <- c(gplots::colorpanel(n = 50, low = '#2f2bad', mid = '#68cded', high = '#f4f4f4'),
                gplots::colorpanel(n = 50, low = '#f4f4f4', mid = "#ff9f19", high = '#d60000'))
  ## Cell type colors
  Clrs <- reactive({
    colors <- c ('blue', 'gold', 'cyan3', 'navy',
                 'forestgreen', 'orange', 'darkolivegreen3',
                 'brown', 'green', 'purple','deepskyblue', 'cyan',
                 'orangered3', 'coral', 'yellow3', "black",
                 'yellow1', 'darkorchid1', 'darksalmon', 'darkseagreen1',
                 'darkslategrey', 'deeppink4', 'green2', 'lemonchiffon1',
                 'lightcyan', 'midnightblue', 'maroon1', 'orange3', 'palegreen',
                 'palevioletred1', 'peru', 'seagreen1', 'red3', 'snow2',
                 'steelblue1', 'turquoise')
    nnodes <- names(chetah$nodetypes[[1]])
    nodes <- paste0("Node", seq_len(length(chetah$nodetypes)))
    M <- max(length(nodes), length(nnodes))
    if (length(nodes) < 24)   {
      gray <- paste0('gray', rev(seq(2, 92, 4)))
    } else if (length(nodes) < 24) {
      gray <- paste0('gray', rev(seq(2, 92, 2)))
    } else {
      gray <- rep('gray70', M)
    }
    if (M > 36) colors <- grDevices::rainbow(M)
    if(!input$colorsplits) {
      names(colors) <- nnodes
      names(gray) <- nodes
    } else {
      names(gray) <- nnodes
      names(colors) <- nodes
    }
    cols <- c(gray, colors)
    cols <- cols[!is.na(names(cols))]
    return(cols)
  })
  ## Reactive event definition
  ReDefineData <- function(event, data, coor) {
    data <- data.frame(data)
    if (!is.null(event)) {
      select <- rownames(coor)[coor[,1] > event[[1]] & coor[,1] < event[[2]] & coor[,2] > event[[3]] & coor[,2] < event[[4]]]
      data <- data[rownames(data) %in% select, , drop = FALSE]
    }
    return(data)
  }
  ## ----------------------------------- Render UI ---------------------------
  ## Extract classification
  classification <- reactive({ Classify(chetah, input$conf_thresh) })
  ## Select a node
  output$nodebutton <- renderUI ({
    req(chetah)
    splitmax <- length(chetah$conf_scores)
    selectInput(inputId = "whichnode",
                label = p("Choose a Node", class = 'opt'),
                choices = seq_len(splitmax))
  })
  ## Select a type
  output$typebutton <- renderUI ({
    req(input$whichnode)
    choices <- names(chetah$nodetypes[[as.numeric(input$whichnode)]])
    selectInput(inputId = "whichtype",
                label = p("Choose a Type", class = 'opt'),
                choices = choices)
  })
  ## Select a type for differential expression
  output$typebutton_DE <- renderUI ({
    req(input$whichnode)
    choices <- names(chetah$nodetypes[[1]])
    selectInput(inputId = "whichtype_DE",
                label = p("Choose a Type", class = 'opt'),
                choices = choices)
  })
  ## Check whether the current type is in the current branch, otherwise stop plotting
  typecorrect <- reactive({
    req(input$whichtype, input$whichnode)
    input$whichtype %in% names(chetah$nodetypes[[as.numeric(input$whichnode)]])
  })
  ## ----------------------------------- Tab 1 ---------------------------
  ## Classification t-SNE
  classTsne <- reactive({
    toplot <- data.frame(classification())
    colnames(toplot) <- 'Cell type'
    PlotTSNE(toplot = toplot, coor = coor, col = Clrs(),
                   pt.size = input$ptsize, return = TRUE, shiny = 'Cell type: ') +
      guides(color=guide_legend(title="Cell types")) +
      guides(color=guide_legend(override.aes = list(size=6)), legend_label = 'Cell type')
  })
  output$classTsne <- plotly::renderPlotly({ plotly::ggplotly(classTsne(), tooltip = 'text') })
  output$dwn_clTsne <- downloadHandler(
    filename = function() { paste0('CHETAH_classification_', input$conf_thresh, '.png') },
    content = function(file) {
      ggsave(file, plot = classTsne(), device = "png", width = 20, height = 12, dpi = 400)
    }
  )
  ## Classification Tree output
  classTree <- reactive({
    if (input$colorsplits) interm <- TRUE else interm <- FALSE
    if(interm) PlotTree(chetah = chetah, col_nodes = Clrs(), no_bgc = T) else {
      PlotTree(chetah = chetah, col = Clrs()) + ggtitle('')
    }
  })
  output$classTree <- renderPlot({ classTree() })
  output$dwn_clTree <- downloadHandler(
    filename = function() { 'CHETAH_tree.png' },
    content = function(file) {
      ggsave(file, plot = classTree(), device = "png", width = 20, height = 12, dpi = 400)
    }
  )
  ## Calculate cell type percentage
  output$typestable <- renderPlot({
    types <- unique(classification())
    types <- c(types[!grepl('Node', types)], sort(types[grepl('Node', types)]))
    perc <- vector()
    for (i in types) perc <- c(perc, round(sum(classification() == i)/length(classification()), 4)*100)
    data <- data.frame(types, perc, stringsAsFactors = TRUE)
    data$types <- factor(types, levels = types)
    ggplot(data, aes_string(x = 'types', y = 'perc', label = 'perc')) +
      geom_bar(stat = 'identity', aes_string(fill = 'types')) +
      theme_classic() +
      scale_fill_manual(values =Clrs()) +
      geom_text(vjust = -0.5, fontface = 'bold')
  })
  ## ----------------------------------- Tab 2 ---------------------------
  ## Plot colored classification tree
  confTree <- reactive({
    ## All black labels
    nodecl <- rep('#000000', length(Clrs()))
    names(nodecl) <- names(Clrs())
    ## Color the branches of the current node
    which <- as.numeric(input$whichnode)
    nodetp <- chetah$nodetypes[[which]]
    lowerns <- (which+1):length(chetah$conf_scores)
    nodes <- unlist(lapply(chetah$nodetypes[lowerns], function (x) {
      if (all(names(x) %in% names(chetah$nodetypes[[which]]))) {
        nodetp[names(x)][1]
      } else 3
    }))
    names(nodes) <- paste0('Node', lowerns)
    nodes <- nodes[nodes !=  3]
    nodecl[c(names(nodetp)[nodetp == 1], names(nodes)[nodes == 1])] <- '#ff0a0a'
    nodecl[c(names(nodetp)[nodetp == 2], names(nodes)[nodes == 2])] <- '#00a1ff'
    ## Plot tree
    PlotTree(chetah, col = nodecl, col_nodes = nodecl, no_bgc = TRUE) + ggtitle('')
  })
  output$confTree <- renderPlot({ confTree() })
  output$dwn_confTree <- downloadHandler(
    filename = function() { paste0('CHETAHtree_node', input$whichnode, '.png') },
    content = function(file) {
      ggsave(file, plot = confTree(), device = "png", width = 20, height = 12, dpi = 400)
    }
  )
  ## Select the max confidence of a cell in the current node
  MaxConf <- reactive ({
    which <- as.numeric(input$whichnode)
    data <- chetah$conf_scores[[which]]
    maxs <- apply(data, 1, max)
    wmax <- apply(data, 1, which.max)
    wmax <- chetah$nodetypes[[which]][colnames(data)][wmax]
    maxs[wmax == 2] <- -maxs[wmax == 2]
    return(maxs)
  })

  ## Select current types incl nodes
  SlctNodes <- reactive ({
    which <- as.numeric(input$whichnode)
    nodes <- unlist(lapply(chetah$nodetypes[(which):length(chetah$conf_scores)], function (x) {
      all(names(x) %in% names(chetah$nodetypes[[which]]))
    }))
    nodes <- paste0('Node', ((which):length(chetah$conf_scores))[nodes])
    nodes <- c(nodes, names(chetah$nodetypes[[which]]))
    return(nodes)
  })

  ## Select the input cell still in the current branch
  SlctCells <- reactive ({
    cellnms <- names(classification())[classification() %in% SlctNodes()]
    return(cellnms)
  })

  ## Confidence t-SNE
  confTsne <- reactive({
    req(SlctCells(), MaxConf(), input$whichnode)
    data <- data.frame(MaxConf()[SlctCells()], classification()[SlctCells()])
    colnames(data) <- c('Confidence', 'cell type')
    coor_prof <- coor[SlctCells(), ]
    coor_prof <- coor_prof[rownames(data), ]
    pl <- PlotTSNE(data, coor_prof,
                   col = grad_col, limits = c(-2,2),
                   pt.size = input$ptsize, return = TRUE, shiny = 'confidence score: ',
                   x_limits = c(min(coor[,1]), max(coor[,1])),
                   y_limits = c(min(coor[,2]), max(coor[,2])),
                   legend_label = 'Confidence')
    pl + labs(color = 'Confidence')
  })
  output$confTsne <- plotly::renderPlotly({ plotly::ggplotly(confTsne(), tooltip = c('text'), source = 'conf') })
  output$dwn_confTsne <- downloadHandler(
    filename = function() { paste0('CHETAH_confidence_node', input$whichnode, '.png') },
    content = function(file) {
      ggsave(file, plot = confTsne(), device = "png", width = 20, height = 12, dpi = 400)
    }
  )

  ## Observe plotly_event, but update to null, when new node is chosen: this doesn't happen automatically
  zoom <- reactive ({ plotly::event_data("plotly_relayout", source = 'conf') })
  rv <- reactiveValues(prof = NULL)
  observeEvent(zoom(), {
    rv$conf <- zoom()
  })
  observeEvent(input$whichnode, {
    rv$conf <- NULL
  })

  ## Confidence heatmap
  confHM <- reactive({
    ## select the data and sort
    which <- as.numeric(input$whichnode)
    data <- chetah$prof_scores[[which]]
    branchtypes <- sort(chetah$nodetypes[[which]])
    data <- data[ ,names(branchtypes)]
    data <- data[names(sort(MaxConf())), ]
    data <- data[rownames(data) %in% SlctCells(), ]
    data1 <- data
    data <- ReDefineData(event = rv$conf, data = data, coor = coor)
    classif <- data.frame('celltypes' = classification()[rownames(data)])

    if (input$sortByType)  data <- data[rownames(classif)[order(classif)], ]
    ann_col <- list('celltypes' = Clrs()[SlctNodes()])
    ## Rescale gradient colors
    rng <- (range(data) + 1)*50
    grad_clrs <- grad_col[rng[1]:rng[2]]

    pheatmap::pheatmap(t(data),
                       color = grad_clrs,
                       border_color = NA,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_colnames = FALSE,
                       scale = "none",
                       annotation_col = classif,
                       annotation_colors = ann_col,
                       height = input$matrixheight)
  })

  output$confHM <- renderPlot({
    if ('RStudioGD' %in% names(dev.list())) dev.off(which = which(names(dev.list()) == 'RStudioGD') + 1)
    print(confHM())
  })

  output$dwn_confHM <- downloadHandler(
    filename = function() { paste0('CHETAH_prof_score_heatmap_node', input$whichnode, '.png') },
    content = function(file) {
      png(file, width = 20, height = 12, res = 400, units = 'in')
      print(confHM())
      dev.off()
    }
  )
  ## ----------------------------------- Tab 3 ---------------------------
  ## Zoom of tSNE must be retained when new data is plotted: this doesn't happen automatically
  zoom2 <- reactive ({ plotly::event_data("plotly_relayout", source = 'prof') })
  rv2 <- reactiveValues(prof = NULL)
  observeEvent(input$whichnode, {
    rv2$prof <- zoom2()
  })
  observeEvent(input$whichtype, {
    rv2$prof <- zoom2()
  })
  observeEvent(zoom2(), {
    rv2$prof <- NULL
  })
  ## Profiel score t-SNE
  profTsne <- reactive({
    req(typecorrect())
    toplot <- data.frame(chetah$prof_scores[[as.numeric(input$whichnode)]][ ,input$whichtype, drop = FALSE],
                         classification())
    colnames(toplot) <- c('Profile score', 'cell type')
    if (!is.null(rv2$prof)) toplot <- ReDefineData(event = zoom2(), data = toplot, coor = coor)
    coor <- coor[rownames(toplot), ]
    pl <- PlotTSNE(toplot, coor,
                   col = grad_col,
                   limits = c(-1,1),
                   pt.size = input$ptsize, return = TRUE, shiny = 'Profile score: ',
                   legend_label = 'Profile score') +
      labs(color = 'Profile score')
    pl
  })

  output$profTsne <- plotly::renderPlotly({ plotly::ggplotly(profTsne(), tooltip = c('text'), source = 'prof') })
  output$dwn_profTsne <- downloadHandler(
    filename = function() { paste0('CHETAH_prof_scores_node', input$whichnode, 'type_', input$whichtype, '.png') },
    content = function(file) {
      png(file, width = 20, height = 12, res = 400, units = 'in')
      print(profTsne())
      dev.off()
    }
  )
  ## Profile score Boxplot
  profBox <- reactive({
    req(typecorrect())
    toplot <- as.data.frame(chetah$prof_scores[[as.numeric(input$whichnode)]][ ,input$whichtype, drop = FALSE])
    clas   <- as.data.frame(chetah$classification)
    event  <- plotly::event_data("plotly_relayout", source = 'prof')
    toplot <- ReDefineData(event = event, data = toplot, coor = coor)
    clas   <- ReDefineData(event = event, data = clas, coor = coor)
    PlotBox(class = clas,
            toplot = toplot,
            col = Clrs())
  })

  output$profBox <- renderPlot({ profBox() })
  output$dwn_profBox <- downloadHandler(
    filename = function() { paste0('CHETAHboxpl_prof_scores_node', input$whichnode, 'type_', input$whichtype, '.png') },
    content = function(file) {
      png(file, width = 20, height = 12, res = 400, units = 'in')
      print(profBox())
      dev.off()
    }
  )

  ## Plot colored classification tree
  celltTree <- reactive({
    req(typecorrect())
    function(){
      ## All black labels
      nodecl <- rep('#000000', length(Clrs()))
      names(nodecl) <- names(Clrs())
      ## Color the current type in the current node
      which <- as.numeric(input$whichnode)
      lowerns <- (which+1):length(chetah$nodetypes)
      nodetp <- chetah$nodetypes[[which]]
      nodes <- unlist(lapply(chetah$nodetypes[lowerns], function (x) {
        if (all(names(x) %in% names(chetah$nodetypes[[which]]))) {
          nodetp[names(x)][1]
        } else 3
      }))
      names(nodes) <- paste0('Node', lowerns)
      nodes <- nodes[nodes !=  3]
      if (input$whichtype %in% names(nodetp)[nodetp == 1]) {
        nodecl[c(names(nodetp)[nodetp == 2], names(nodes)[nodes == 2])] <- '#ff0a0a'
        nodecl[input$whichtype] <- '#00a1ff'
      } else {
        nodecl[c(names(nodetp)[nodetp == 1], names(nodes)[nodes == 1])] <- '#00a1ff'
        nodecl[input$whichtype] <- '#ff0a0a'
      }
      ## Plot tree
      PlotTree(chetah, col = nodecl, col_nodes = nodecl, no_bgc = TRUE) + ggtitle('')
    }
  })
  output$profTree <- renderPlot({ celltTree()() })
  ## ----------------------------------- Tab4  ---------------------------
  ## Gene heatmap
  exprHM <- reactive({
    req(typecorrect())
    ## select the data and sort
    which <- as.numeric(input$whichnode)
    ## Cells
    branchtypes <- sort(chetah$nodetypes[[which]])
    cells <- classification()[classification() %in% names(branchtypes)]
    annotation_col <- data.frame('celltypes' = cells)
    cells <- cells[order(match(cells, names(branchtypes)))]
    ## genes
    if (input$largediff) {
      lh_genes <- chetah$genes[[which]][[input$whichtype]]
      data <- as.matrix(counts[names(lh_genes), names(cells)])
      r_branch <- apply(data[ ,names(cells)[cells %in% names(branchtypes)[branchtypes == 1]]], 1, mean)
      l_branch <- apply(data[ ,names(cells)[cells %in% names(branchtypes)[branchtypes == 2]]], 1, mean)
      lh_genes <- lh_genes[names(sort(abs(r_branch - l_branch), decreasing = TRUE))[seq_len(input$n_genes)]]
      lh_genes <- sort(lh_genes)
    } else {
      lh_genes <- sort(chetah$genes[[which]][[input$whichtype]][seq_len(input$n_genes)])
    }
    lh_genes[lh_genes == 0] <- 'low'
    lh_genes[lh_genes == 1] <- 'high'
    annotation_row <- data.frame('High_Low' = lh_genes)
    genes <- names(lh_genes)
    data <- as.matrix(counts[genes, names(cells)])

    hl_clrs <- c('black', 'gray80')
    names(hl_clrs) <- c('low', 'high')
    ann_col <- list('celltypes' = Clrs()[names(branchtypes)],
                    'High_Low' = hl_clrs)
    if (input$scaling) scale <- 'row' else scale <- 'none'
    ## Plot
    pheatmap::pheatmap(data,
                       color = grad_col,
                       border_color = NA,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_colnames = FALSE,
                       scale = scale,
                       annotation_col = annotation_col,
                       annotation_row = annotation_row,
                       annotation_colors = ann_col,
                       fontsize_row = input$lettersize)
  })
  output$exprHM <- renderPlot({
    if ('RStudioGD' %in% names(dev.list())) dev.off(which = which(names(dev.list()) == 'RStudioGD') + 1)
    exprHM()
  })

  output$dwn_exprHM <- downloadHandler(
    filename = function() { paste0('CHETAHgenes_node', input$whichnode, 'type_', input$whichtype, '.png') },
    content = function(file) {
      png(file, width = 20, height = 12, res = 400, units = 'in')
      print(exprHM())
      dev.off()
    }
  )
  ## Plot classification tree
  output$genesTree <- renderPlot({ celltTree()() })
  ## ----------------------------------- Tab 5  ---------------------------
  fld_ch <- reactive({
    type <- input$whichtype_DE
    otherbranch <- names(chetah$nodetypes[[1]])[names(chetah$nodetypes[[1]]) != type]
    apply(counts[ ,names(classification())[classification() %in% otherbranch], drop = FALSE], 1, mean) -
                apply(counts[ ,names(classification())[classification() == type], drop = FALSE], 1, mean)
  })
  ## Gene heatmap
  diff_exp_HM <- reactive({
    req(input$whichtype_DE)
    ## Select types and cells
    type <- input$whichtype_DE
    otherbranch <- names(chetah$nodetypes[[1]])[names(chetah$nodetypes[[1]]) != type]
    branchtypes <- c(type, otherbranch)
    cells <- classification()[classification() %in% branchtypes]
    cells <- cells[order(match(cells, branchtypes))]
    ## Select genes + sort
    fld_genes <-  names(sort(abs(fld_ch()), decreasing = TRUE)[1:input$n_genes])
    lh_genes <- rep('high', input$n_genes)
    names(lh_genes) <- fld_genes
    lh_genes[fld_ch()[fld_genes] < 0] <- 'low'
    lh_genes <- sort(lh_genes)
    ## Select data
    data <- counts[names(lh_genes), names(cells)]
    ## Annotation
    annotation_col <- data.frame('celltypes' = cells)
    annotation_row <- data.frame('High_Low' = lh_genes)
    hl_clrs <- c('black', 'gray80')
    names(hl_clrs) <- c('low', 'high')
    ann_col <- list('celltypes' = Clrs()[branchtypes],
                    'High_Low' = hl_clrs)
    if (input$scaling) scale <- 'row' else scale <- 'none'
    ## Plot
    pheatmap::pheatmap(data,
                       color = grad_col,
                       border_color = NA,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_colnames = FALSE,
                       scale = scale,
                       annotation_col = annotation_col,
                       annotation_row = annotation_row,
                       annotation_colors = ann_col,
                       fontsize_row = input$lettersize)
  })
  output$diff_exp_HM <- renderPlot({
    if ('RStudioGD' %in% names(dev.list())) dev.off(which = which(names(dev.list()) == 'RStudioGD') + 1)
    diff_exp_HM()
  })
  output$dwn_diff_exp_HM <- downloadHandler(
    filename = function() { paste0('diff_expr_genes_type', input$whichtype, '.png') },
    content = function(file) {
      png(file, width = 20, height = 12, res = 400, units = 'in')
      print(diff_exp_HM())
      dev.off()
    }
  )
}

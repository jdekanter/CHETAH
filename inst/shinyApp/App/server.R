server <- function(input, output, session) {
  ## ----------------------------------- Startup notification ------------------
  observe({
    showNotification("Welcome to the CHETAH shiny application.
Here, you can interact with your data and the CHETAH output.
Most plots have a help-message that can be viewed by hovering over the title of the plot.
Download buttons are present below each plot.
(These functions work better in a web-browser: click 'Open in Browser' on top of the Rstudio-window.)
Each t-SNE plot is interactive, so you can zoom in and out. Other plots will adapt
to only show information of the cells in the currently selected t-SNE window.
See the 'Info' tab for info about the method.", duration = 60)
  })

  ## Make short links (variables) to the chetah data [links without adjustment don't take up space]
  coor <- SingleCellExperiment::reducedDim(CHETAH:::ch_env$chetah, CHETAH:::ch_env$redD)
  counts <- SummarizedExperiment::assay(CHETAH:::ch_env$chetah, CHETAH:::ch_env$input_c)
  conf_scores <- CHETAH:::ch_env$chetah@int_colData$CHETAH$conf_scores
  prof_scores <- CHETAH:::ch_env$chetah@int_colData$CHETAH$prof_scores
  meta_data <- CHETAH:::ch_env$chetah@int_metadata$CHETAH
  chetah <- CHETAH:::ch_env$chetah
  redD <- CHETAH:::ch_env$redD
  input_c <- CHETAH:::ch_env$input_c

  ## ----------------------------------- Select colors ------------------
  ## Confidence t-SNE color gradient
  grad_col <- c(gplots::colorpanel(n = 50, low = '#2f2bad', mid = '#68cded', high = '#f4f4f4'),
                gplots::colorpanel(n = 50, low = '#f4f4f4', mid = "#ff9f19", high = '#d60000'))
  ## Cell type colors
  Clrs <- reactive({
    req(chetah)
    colors <- c ('blue', 'gold', 'cyan3', 'navy',
                 'forestgreen', 'orange', 'darkolivegreen3',
                 'brown', 'green', 'purple','deepskyblue', 'cyan',
                 'orangered3', 'coral', 'yellow3', "black",
                 'yellow1', 'darkorchid1', 'darksalmon', 'darkseagreen1',
                 'darkslategrey', 'deeppink4', 'green2', 'lemonchiffon1',
                 'lightcyan', 'midnightblue', 'maroon1', 'orange3', 'palegreen',
                 'palevioletred1', 'peru', 'seagreen1', 'red3', 'snow2',
                 'steelblue1', 'turquoise')
    nnodes <- names(meta_data$nodetypes[[1]])
    nodes <- c("Unassigned", paste0("Node", seq_len(length(meta_data$nodetypes) - 1)))

    if(input$colornodes) {
      n <- nodes
      nodes <- nnodes
      nnodes <- n
      rm(n)
    }
    if (length(nodes) < 24)   {
      gray <- paste0('gray', rev(seq(2, 92, 4)))
    } else if (length(nodes) < 24) {
      gray <- paste0('gray', rev(seq(2, 92, 2)))
    } else {
      gray <- rep('gray70', length(nodes))
    }
    if (length(nnodes) > 36) colors <- grDevices::rainbow(M)
    names(colors) <- nnodes
    names(gray) <- nodes
    cols <- c(gray, colors)
    cols <- cols[!is.na(names(cols))]
    return(cols)
  })
  ## Reactive event definition
  ReDefineData <- function(event, data, coor) {
    data <- data.frame(data)
    if (!is.null(event)) {
      if(names(event)[1] == "xaxis.range[0]") {
        select <- rownames(coor)[coor[,1] > event[["xaxis.range[0]"]] &
                                   coor[,1] < event[["xaxis.range[1]"]] &
                                   coor[,2] > event[["yaxis.range[0]"]] &
                                   coor[,2] < event[["yaxis.range[1]"]]]
        data <- data[rownames(data) %in% select, , drop = FALSE]
      }
    }
    return(data)
  }
  ## ----------------------------------- Render UI ---------------------------
  ## Extract classification
  classification <- reactive({
    req(input$conf_thresh, chetah)
    Classify(chetah, input$conf_thresh, return_clas = TRUE)
  })
  ## Select a node
  output$nodebutton <- renderUI ({
    req(chetah)
    splitmax <- length(conf_scores)
    choices <- seq_len(splitmax)
    names(choices) <- c(0, seq_len(splitmax - 1))
    selectInput(inputId = "whichnode",
                label = p("Choose a Node", class = 'opt'),
                choices = choices)
  })
  ## Select a type
  output$typebutton <- renderUI ({
    req(input$whichnode)
    choices <- names(meta_data$nodetypes[[as.numeric(input$whichnode)]])
    selectInput(inputId = "whichtype",
                label = p("Choose a Type", class = 'opt'),
                choices = choices)
  })
  ## Select a type for differential expression
  output$typebutton_DE <- renderUI ({
    req(input$whichnode)
    choices <- names(meta_data$nodetypes[[1]])
    selectInput(inputId = "whichtype_DE",
                label = p("Choose a Type", class = 'opt'),
                choices = choices)
  })
  ## Check whether the current type is in the current branch, otherwise stop plotting
  typecorrect <- reactive({
    req(input$whichtype, input$whichnode)
    input$whichtype %in% names(meta_data$nodetypes[[as.numeric(input$whichnode)]])
  })
  ## Select a gene for Tab 5
  output$geneselection <- renderUI ({
    req(counts)
    choices <- rownames(counts)
    selectInput(inputId = "whichgene",
                label = p("Choose a Gene", class = 'opt'),
                choices = choices)
  })
  ## ----------------------------------- Tab 1 ---------------------------
  ## Classification t-SNE
  classTsne <- reactive({
    req(classification(), coor)
    toplot <- classification()
    toplot[toplot == "Unassigned"] <- "Unassigned (Node0)"
    u_toplot <- unique(toplot)
    toplot <- data.frame(factor(toplot, levels = c(u_toplot[grepl("Unassigned", u_toplot)],
                                        sort(u_toplot[grepl("Node", u_toplot) & !grepl("\\(Node0)", u_toplot)]),
                                        sort(u_toplot[!grepl("Node|Unassigned", u_toplot)]))))
    colnames(toplot) <- 'Cell type'
    clrs <- Clrs()
    names(clrs)[names(clrs) == "Unassigned"] <- "Unassigned (Node0)"
    PlotTSNE(toplot = toplot, input = chetah, redD = redD, col = clrs,
             pt.size = input$ptsize, return = TRUE, shiny = 'Cell type: ') +
      guides(color=guide_legend(title="Cell types")) +
      guides(color=guide_legend(override.aes = list(size=6)), legend_label = 'Cell type')
  })
  output$classTsne <- plotly::renderPlotly({
    classTsne_plot <- plotly::ggplotly(classTsne(), tooltip = 'text')
    plotly::config(classTsne_plot, modeBarButtonsToRemove = c('toggleSpikelines','lasso2d', 'select2d', 'hoverCompareCartesian'), collaborate = FALSE)
  })
  output$dwn_clTsne <- downloadHandler(
    filename = function() { paste0('CHETAH_classification_', input$conf_thresh, '.png') },
    content = function(file) {
      ggsave(file, plot = classTsne(), device = "png", width = 20, height = 12, dpi = 400)
    }
  )
  ## Classification Tree output
  classTree <- reactive({
    req(Clrs())
    if (input$colornodes) interm <- TRUE else interm <- FALSE
    if(interm) PlotTree(input = chetah, col_nodes = Clrs(), no_bgc = TRUE) else {
      PlotTree(input = chetah, col = Clrs()) + ggtitle('')
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
  typestable <- reactive({
    req(classification())
    types <- unique(classification())
    types <- c(sort(types[!grepl('Node|Unassigned', types)]), types[grepl('Unassigned', types)] , sort(types[grepl('Node', types)]))
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
  output$typestable <- renderPlot({ typestable() })
  output$dwn_ttable <- downloadHandler(
    filename = function() { 'CHETAH_percentages.png' },
    content = function(file) {
      ggsave(file, plot = typestable(), device = "png", width = 20, height = 12, dpi = 400)
    }
  )
  ## ----------------------------------- Tab 2 ---------------------------
  ## Plot colored classification tree
  confTree <- reactive({
    req(Clrs(), input$whichnode)
    ## All black labels
    nodecl <- rep('#000000', length(Clrs()))
    names(nodecl) <- names(Clrs())
    ## Color the branches of the current node
    which <- as.numeric(input$whichnode)
    nodetp <- meta_data$nodetypes[[which]]
    lowerns <- (which+1):length(conf_scores)
    nodes <- unlist(lapply(meta_data$nodetypes[lowerns], function (x) {
      if (all(names(x) %in% names(meta_data$nodetypes[[which]]))) {
        nodetp[names(x)][1]
      } else 3
    }))
    names(nodes) <- paste0('Node', lowerns - 1) ## The names are shifted 1 place
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
    req(input$whichnode)
    which <- as.numeric(input$whichnode)
    data <- conf_scores[[which]]
    maxs <- apply(as.matrix(data), 1, max)
    wmax <- apply(as.matrix(data), 1, which.max)
    wmax <- meta_data$nodetypes[[which]][colnames(data)][wmax]
    maxs[wmax == 2] <- -maxs[wmax == 2]
    return(maxs)
  })

  ## Select current types incl nodes
  SlctNodes <- reactive ({
    req(input$whichnode)
    which <- as.numeric(input$whichnode)
    nodes <- unlist(lapply(meta_data$nodetypes[(which):length(conf_scores)], function (x) {
      all(names(x) %in% names(meta_data$nodetypes[[which]]))
    }))
    nodes <- paste0('Node', ((which - 1):(length(conf_scores) - 1))[nodes])
    nodes[grepl("Node0", nodes)] <- "Unassigned"
    nodes <- c(nodes, names(meta_data$nodetypes[[which]]))
    return(nodes)
  })

  ## Select the input cell still in the current branch
  SlctCells <- reactive ({
    req(SlctNodes, classification())
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
    pl <- PlotTSNE(toplot = data, input = chetah, redD = redD,
                   col = grad_col, limits = c(-2,2),
                   pt.size = input$ptsize, return = TRUE, shiny = 'confidence score: ',
                   x_limits = c(min(coor[,1]), max(coor[,1])),
                   y_limits = c(min(coor[,2]), max(coor[,2])),
                   legend_label = 'Confidence')
    pl + labs(color = 'Confidence')
  })
  output$confTsne <- plotly::renderPlotly({
    confTsne_plot <- plotly::ggplotly(confTsne(), tooltip = c('text'), source = 'conf')
    plotly::config(confTsne_plot, modeBarButtonsToRemove = c('toggleSpikelines','lasso2d', 'select2d', 'hoverCompareCartesian'), collaborate = FALSE)
  })
  output$dwn_confTsne <- downloadHandler(
    filename = function() { paste0('CHETAH_confidence_node', input$whichnode, '.png') },
    content = function(file) {
      ggsave(filename = file, plot = confTsne(), device = "png", width = 20, height = 12, dpi = 400)
    }
  )

  ## Observe plotly_event, but update to null, when new node is chosen: this doesn't happen automatically
  zoom <- reactive ({ plotly::event_data("plotly_relayout", source = 'conf') })
  rv <- reactiveValues(conf = NULL)
  observeEvent(zoom(), {
    rv$conf <- zoom()
  })
  observeEvent(input$whichnode, {
    rv$conf <- NULL
  })

  ## Confidence heatmap
  confHM <- reactive({
    req(SlctCells(), MaxConf(), input$whichnode, rv)
    ## select the data and sort
    which <- as.numeric(input$whichnode)
    data <- prof_scores[[which]]
    branchtypes <- sort(meta_data$nodetypes[[which]])
    data <- data[ ,names(branchtypes)]
    data <- data[names(sort(MaxConf())), ]
    data <- data[rownames(data) %in% SlctCells(), ]
    data <- ReDefineData(event = rv$conf, data = data, coor = coor)
    classif <- data.frame('celltypes' = classification()[rownames(data)])
    if (nrow(data) > 0) {
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
    } else {
      NULL
    }
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
    req(typecorrect(), chetah, rv2)
    toplot <- data.frame(prof_scores[[as.numeric(input$whichnode)]][ ,input$whichtype, drop = FALSE],
                         classification())
    colnames(toplot) <- c('Profile score', 'cell type')
    if (!is.null(rv2$prof)) toplot <- ReDefineData(event = zoom2(), data = toplot, coor = coor)
    coor <- coor[rownames(toplot), ]
    pl <- PlotTSNE(toplot, input = chetah, redD = redD,
                   col = grad_col,
                   limits = c(-1,1),
                   pt.size = input$ptsize, return = TRUE, shiny = 'Profile score: ',
                   legend_label = 'Profile score') +
      labs(color = 'Profile score')
    pl
  })
  output$profTsne <- plotly::renderPlotly({
    profTsne_plot <- plotly::ggplotly(profTsne(), tooltip = c('text'), source = 'prof')
    plotly::config(profTsne_plot, modeBarButtonsToRemove = c('toggleSpikelines','lasso2d', 'select2d', 'hoverCompareCartesian'), collaborate = FALSE)
  })
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
    req(typecorrect(), chetah)
    toplot <- as.data.frame(prof_scores[[as.numeric(input$whichnode)]][ ,input$whichtype, drop = FALSE])
    clas   <- as.data.frame(classification())
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
      lowerns <- (which+1):length(meta_data$nodetypes)
      nodetp <- meta_data$nodetypes[[which]]
      nodes <- unlist(lapply(meta_data$nodetypes[lowerns], function (x) {
        if (all(names(x) %in% names(meta_data$nodetypes[[which]]))) {
          nodetp[names(x)][1]
        } else 3
      }))
      names(nodes) <- paste0('Node', lowerns - 1)
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
  ## Select cells
  HMcells <- reactive({
    req(typecorrect(), classification(), SlctNodes())
    which <- as.numeric(input$whichnode)
    ## Select cells
    branchtypes <- names(sort(meta_data$nodetypes[[which]]))
    if (input$inclnodes) {
      branchtypes <- c(branchtypes, SlctNodes()[grepl("Node|Unassigned", SlctNodes())])
    }
    cells <- classification()[classification() %in% branchtypes]
    cells <- cells[order(match(cells, branchtypes))]
    validate(
      need(cells > 0, "
           No cells to display in this Node. Lowering the 'Confidence Threshold' might work.
           Please check the 'Classification' tab to see the classification with the current 'Confidence Threshold'")
      )
    return(cells)
  })

  ## Select genes
  HMgenes <- reactive({
    req(HMcells())
    which <- as.numeric(input$whichnode)
    ## Select Genes
    if (input$largediff) {
      finaltypes <- sort(meta_data$nodetypes[[which]])
      lh_genes <- meta_data$genes[[which]][[input$whichtype]]
      data <- as.matrix(counts[names(lh_genes), names(HMcells())])
      r_branch <- apply(as.matrix(data[ ,names(HMcells())[HMcells() %in% names(finaltypes)[finaltypes == 1]]]), 1, mean)
      l_branch <- apply(as.matrix(data[ ,names(HMcells())[HMcells() %in% names(finaltypes)[finaltypes == 2]]]), 1, mean)
      lh_genes <- lh_genes[names(sort(abs(r_branch - l_branch), decreasing = TRUE))[seq_len(input$n_genes)]]
      lh_genes <- sort(lh_genes)
      ## If no cells in one of the two branches: Don't sort
      if (length(lh_genes) == 0) lh_genes <- sort(meta_data$genes[[which]][[input$whichtype]][seq_len(input$n_genes)])
    } else {
      lh_genes <- sort(meta_data$genes[[which]][[input$whichtype]][seq_len(input$n_genes)])
    }
    lh_genes[lh_genes == 0] <- 'lowly'
    lh_genes[lh_genes == 1] <- 'highly'
    return(lh_genes)
  })

  ## Gene heatmap
  exprHM <- reactive({
    req(HMgenes())
    which <- as.numeric(input$whichnode)
    ## Select cells
    branchtypes <- names(sort(meta_data$nodetypes[[which]]))
    if (input$inclnodes) {
      branchtypes <- c(branchtypes, SlctNodes()[grepl("Node|Unassigned", SlctNodes())])
    }
    genes <- names(HMgenes())
    data <- as.matrix(counts[genes, names(HMcells())])

    hl_clrs <- c('black', 'gray80')
    names(hl_clrs) <- c('lowly', 'highly')
    annotation_col <- data.frame('celltypes' = HMcells())
    annotation_row <- data.frame('Expressed_in_this_type' = HMgenes())
    ann_col <- list('celltypes' = Clrs()[branchtypes],
                    'Expressed_in_this_type' = hl_clrs)
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
                       fontsize_row = input$lettersize,
                       annotation_names_row = FALSE)
  })
  output$exprHM <- renderPlot({
    if ('RStudioGD' %in% names(dev.list())) dev.off(which = which(names(dev.list()) == 'RStudioGD') + 1)
    exprHM()
  })

  output$dwn_exprHM <- downloadHandler(
    filename = function() { paste0('CHETAH_HM_node', input$whichnode, '_', input$whichtype, '.png') },
    content = function(file) {
      png(file, width = 20, height = 12, res = 400, units = 'in')
      print(exprHM())
      dev.off()
    }
  )
  ## Download the genes!
  output$dwn_genes <- downloadHandler(
    filename = function() { paste0('CHETAHgenes_node', input$whichnode, '_', input$whichtype, '.txt') },
    content = function(con) {
      write.table(matrix(names(HMgenes()), ncol = 1), file = con, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  )
  ## Plot classification tree
  output$genesTree <- renderPlot({ celltTree()() })
  ## ----------------------------------- Tab5  ---------------------------
  geneExpr <- reactive({
    req(input$whichgene, classification())
    genes <- counts[input$whichgene, ]
    #if(!is.null(log)) genes <- log2(genes + 1)
    genes <- cbind.data.frame(genes, as.data.frame(classification()))
    #if(!is.null(sub)) genes <- genes[genes$type %in% sub, ]
    colnames(genes) <- c('gene', 'type')
    suppressMessages(genes <- reshape2::melt(genes))

    ggplot(genes, aes(x=type, y = value)) +
      geom_jitter(aes(color = type), size = 2) +
      geom_boxplot(aes(color = type), fill = rgb(1,1,1,0.5), outlier.colour = "NA")  +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
            axis.title.x = element_blank()) +
      guides(color=FALSE) +
      scale_color_manual(values = Clrs())
  })
  output$geneExpr <- renderPlot({ geneExpr() })
  output$dwn_geneExpr <- downloadHandler(
    filename = function() { paste0('CHETAH_gene_exp_', input$whichgene, '.png') },
    content = function(file) {
      png(file, width = 20, height = 12, res = 400, units = 'in')
      print(geneExpr())
      dev.off()
    }
  )
}

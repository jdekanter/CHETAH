# Define UI ----
ui <- fluidPage(
  tags$head(

    ## CSS Style
    tags$style(HTML("
h1 {
  font-family: 'Raleway Black', 'Arial Black', sans-serif;
  color: #ffffff;
}
.options {
  background: #000000;
  min-height: 100%;
}
.opt {
  font-family: 'Raleway ExtraBold', 'Arial Black', sans-serif;
  font-weight: 'bold';
  color: #fcb22a;
}
.value {
  font-family: 'Raleway', 'Arial', sans-serif;
}
h3 {
  font-family: 'Raleway ExtraBold', 'Arial Black', sans-serif;
}
.tabbable > .nav > li > a                  {background-color: #fef2b9;  color:black; font-family: 'Raleway ExtraBold', 'Arial Black', sans-serif;}
.tabbable > .nav > li[class=active]    > a {background-color: #fcb22a;  color:black; font-family: 'Raleway ExtraBold', 'Arial Black', sans-serif;}
                    "))
  ),
  fluidRow(
    #fixedPanel(
      img(src="Shiny_bar.png", align = "left", height = 106, width = 1188)
    #, style="z-index:1;")
  ),
  #
  #
  fluidRow(
    ## OPTION PANEL 1
    column(2,
           fixedPanel(class = 'options',
                      width = "16%",
                      left = '0px',
                      top = '106px',
                      h1("OPTIONS"),
                      ## Color the intermediate types?
                      checkboxInput(inputId = "colornodes",
                                    label = p('Color the intermediate types', class = 'opt'),
                                    value = FALSE),
                      ## Conf threshhold for which cells to classify
                      numericInput("conf_thresh", label = p("Confidence threshold", class = 'opt'),
                                   value = 0.1, step = 0.01),

                      conditionalPanel(
                        condition = "input.MP != 'Genes used by CHETAH' && input.MP != 'Differentially expressed genes' && input.MP != 'Expression per gene'",
                        ## t-SNE point size
                        numericInput("ptsize", label = p("Point Size", class = 'opt'),
                                     value = 1.5, step = 0.5)
                      ),
                      conditionalPanel(
                        condition = "input.MP != 'Classification' && input.MP != 'Differentially expressed genes' && input.MP != 'Expression per gene'",
                        ## Which node to display 2
                        uiOutput(outputId = "nodebutton")
                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Profile scores' || input.MP == 'Genes used by CHETAH'",
                        ## Which type to display 3
                        uiOutput(outputId = "typebutton")
                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Genes used by CHETAH' || input.MP == 'Differentially expressed genes'",
                        ## Scale the matrix 3
                        checkboxInput(inputId = "scaling",
                                      label = p('Scale Matrix', class = 'opt'),
                                      value = FALSE),
                        ## Number of genes 3
                        numericInput(inputId = "n_genes",
                                     label = p('# of genes', class = 'opt'),
                                     value = 200, step = 1, min = 2, max = 200),
                        ## Height of the matrix 3
                        numericInput(inputId = "lettersize",
                                     label = p('Matrix label size', class = 'opt'),
                                     value = 10, step = 0.5, min = 1, max = 20)
                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Differentially expressed genes'",
                        ## Which type to display in 4
                        uiOutput(outputId = "typebutton_DE")

                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Genes used by CHETAH'",
                        ## Height of the matrix 3
                        checkboxInput(inputId = "largediff",
                                      label = p('Genes with max. difference in the INPUT', class = 'opt'),
                                      value = FALSE)
                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Expression per gene'",
                        ## Height of the matrix 3
                        uiOutput(outputId = 'geneselection')
                      )
           )
    ),
    ## MAIN PANEL
    column(10,
           fluidRow(
             tabsetPanel(id = "MP",
                         ### Tab 1
                         tabPanel("Classification",
                                  fluidRow(
                                    column(12,
                                           ## Classification t-SNE
                                           h3('CHETAH classification'),
                                           plotly::plotlyOutput(outputId = 'classTsne', height = 800),
                                           downloadButton('dwn_clTsne', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Classification percentages
                                           h3("Cell type percentages"),
                                           plotOutput('typestable'),
                                           downloadButton('dwn_ttable', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Classification tree
                                           h3("Classification tree"),
                                           plotOutput(outputId = "classTree", height = 600),
                                           downloadButton('dwn_clTree', 'Download Plot')
                                    )
                                  )
                         ),
                         ### Tab 2
                         tabPanel("Confidence scores",
                                  fluidRow(
                                    column(12,
                                           ## Confidence plot
                                           h3('Confidence scores per nodes'),
                                           plotly::plotlyOutput(outputId = "confTsne", height = 600),
                                           downloadButton('dwn_confTsne', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(4,
                                           h3('Current branches'),
                                           ## Colored cl. tree
                                           plotOutput(outputId = "confTree", height = 600),
                                           downloadButton('dwn_confTree', 'Download Plot')
                                    ),
                                    column(8,
                                           h3('Profile score heatmap'),
                                           ## Heatmap
                                           plotOutput(outputId = "confHM", height = 600),
                                           checkboxInput(inputId = "sortByType",
                                                         label = p('Sort by Type', class = 'opt'),
                                                         value = FALSE),
                                           downloadButton('dwn_confHM', 'Download Plot')
                                    )
                                  )
                         ),
                         ### Tab 3
                         tabPanel("Profile scores",
                                  fluidRow(
                                    column(12,
                                           ## Profile score plots
                                           h3('Profiles scores per nodes per type'),
                                           plotly::plotlyOutput(outputId = "profTsne", height = 600),
                                           downloadButton('dwn_profTsne', 'Download Plot')
                                           )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Profile scores in boxplots
                                           h3('In a boxplot'),
                                           plotOutput(outputId = "profBox", height = 600),
                                           downloadButton('dwn_profBox', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Tree
                                           h3('Tree'),
                                           plotOutput(outputId = "profTree", height = 600)
                                    )
                                  )
                         ),
                         ### Tab 4
                         tabPanel("Genes used by CHETAH",
                                  fluidRow(
                                    column(12,
                                           ## Gene heatmap
                                           h3("Selected genes per cell type"),
                                           plotOutput(outputId = "exprHM", height = 800),
                                           downloadButton('dwn_exprHM', 'Download Plot'),
                                           h3("Tree"),
                                           plotOutput(outputId = 'genesTree', height = 600)
                                    )
                                  )
                         ),
                         ### Tab 5
                         tabPanel("Expression per gene",
                                  fluidRow(
                                    column(12,
                                           ## Gene heatmap
                                           h3("Gene expression per cell type"),
                                           plotOutput(outputId = "geneExpr", height = 400),
                                           downloadButton('dwn_geneExpr', 'Download Plot')
                                    )
                                  )
                         ),
                         ### Tab0
                         ### Tab 6
                         tabPanel("Info",
                                  tags$div(
                                    tags$p("Welcome to the CHETAH shiny application."),
                                    tags$p("Here, you can interact with your data and the CHETAH output."),
                                    tags$p("On this page, you will find a short explanation of every tab panel on this browser page,"),
                                    tags$p("the plots on every page and the way to change the plots."),
                                    tags$p("A button below each plot gives the option to download the current plot."),
                                    tags$p('Each t-SNE plot is interactive, so you can zoom in and out. Other plots will adapt to only show information of the cells in the currently selected t-SNE window'),
                                    tags$h3('General'),
                                    tags$p('CHETAH creates a classifcation tree, by hierarchically clustering the reference data.'),
                                    tags$p('CHETAh calculates two types of scores for each input cell j in each node of the classification tree.'),
                                    tags$p('First, for each cell j, a Profile Score is calculated for each reference cell type in this node'),
                                    tags$p('This Profile Score is a -1:1 value that represents the similarity of j to the reference profile.'),
                                    tags$p('Based of these Profile Scores, cell j is assigned to the right or the left branch of the node.'),
                                    tags$p('For this assignment, one Confidence Score is calculated. This score (values of 0:2) represents the confidence of the assignment'),
                                    tags$p('If the assignment of j has a confidence below the Confidence Threshold, the cell is not assigned to a branch, but for j, the classification will stop in that node.'),
                                    tags$p("A classifcation to a non-leaf node are called an 'intermediate type'."),
                                    tags$p("Classifcations to a leaf node of the tree is called a 'final type'."), tags$br(),
                                    tags$p("Recap: for each cell j, in each node, a Profile Score is calculated for each reference profile, the cell is assigned to a branch, and for this assignment, one Confidence Score is calculated."),
                                    tags$h3('Classification'),
                                    tags$p("This is the page where you can view the CHETAH classification."),
                                    tags$p("The first plot shows the classification on a t-SNE map with the coordinates that you provided."),
                                    tags$p('Each plot is colored according to the CHETAH infered cell type'),
                                    tags$p("The second plot shows the percentages of cell types in you sample in a barplot."),
                                    tags$p("The third plot shows the classification tree that CHETAH constructed. The nodes are colored with the same colors as used in the plots above."), tags$br(),
                                    tags$p("To color the intermediate types (/nodes) instead of the final types, click the button in the left menu panel."),
                                    tags$p("To change the strictness of the classifcation, change the Confidence Threshold. The Confidence score has a value between 0 and 2."),
                                    tags$p("0.1 is the default threshold. A threshold of 0 will classify all cells to a final type (leaf node of the tree) and 1 is a very stringent cut-off,"),
                                    tags$p("in which case only the cells with very high confidence will be classified."),
                                    tags$h3('Confidence scores'),
                                    tags$p('Again, the first plot is a t-SNE, but now, the cells are colored by the confidence score of the currently selected node.'),
                                    tags$p('If a cell is assigned to the left branch, the scores are negative (for plotting purposes) and the colors shades of blue. For the right branch shades of red are used.'),
                                    tags$p('The bottom left plot shows the classification tree. Blue cells in the t-SNE plot above are assigned to the blue branches, the red cells are assigned to the red branches.'),
                                    tags$p('The plot in the right corner shows a heatmap of the profile scores for each individual reference profile in the current branch.'),
                                    tags$p('For an input cell j, the more strongly the colors/scores differ between types of the two branches, the higher the confidence score for that cell is'),
                                    tags$br(),
                                    tags$p("The branch that is shown can be changed by the 'Choose Node' bar"),
                                    tags$h3('Profile scores'),
                                    tags$p("All plots in this tab show the profile scores for one reference cell type in one node."),
                                    tags$p("This score represents the chance of a cell being of that type, rather than the types in the other branch."),
                                    tags$p("The first plot shows a t-SNE colored by the profile score currently selected reference profile."),
                                    tags$p("The second plot shows the same values in box plots, grouped by the CHETAH cell type."),
                                    tags$p("The third plot shows the tree, but now the selected reference profile is colored in one color and the types of the other branch are colored in another color. These are the types with which the current type is compared. Is a Profile Score high for type X (with Y and Z in the other branch), than this means that this cell is highly more likely to be of type X than of the types in the other branch."),
                                    tags$br(),
                                    tags$p('In all plots in this tab, only one profile score is depicted. In the left menu, both the node and the reference profile for which the reference profiel should be depicted can be selected.'),
                                    tags$h3('Genes used by CHETAH'),
                                    tags$p('In this panel, the expression in the input data of the genes that were used by CHETAH for the classification can be viewed in a heatmap.'),
                                    tags$p("CHETAH selects a separate group of genes for each reference type, in each node. The cell types to which the current type is compared, can be viewed in the tree below the heatmap. This tree is similar to that on the 'Profile Scores' tab."),
                                    tags$p("The node and reference cell type for which to show the genes can be selected.
                                    Also, the matrix can be scaled per row (gene), by clicking 'Scale Matrix'.
                                    Moreover, the number of genes to show can be adjusted by '# of genes'.
                                    Automatically, the genes that had the highest difference in the reference will be selected.
                                    If the genes that have the highest difference in the input data should be selected, check the 'Genes with max. difference in input' box")
                                  )
                            )
             )
           )
    )
  )
)

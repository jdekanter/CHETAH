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
  tags$head(tags$style(
    "body { word-wrap: break-word; }"
  )),
  fluidRow(
      img(src="Shiny_bar.png", align = "left", height = "80vh", width = "820vh")
  ),
  #
  #
  fluidRow(
    ## OPTION PANEL 1
    column(2,
           left = '0px',
           fluidRow(style = "background:#000000;",
                      h1("OPTIONS"),
                      ## Color the intermediate types?
                      checkboxInput(inputId = "colornodes",
                                    label = p('Color the intermediate types', class = 'opt'),
                                    value = FALSE),
                      ## Conf threshhold for which cells to classify
                      numericInput("conf_thresh", label = p("Confidence threshold", class = 'opt'),
                                   value = 0.1, step = 0.01),

                      conditionalPanel(
                        condition = "input.MP != 'Genes used by CHETAH' && input.MP != 'Expression per gene'",
                        ## t-SNE point size
                        numericInput("ptsize", label = p("Point Size", class = 'opt'),
                                     value = 1.5, step = 0.5)
                      ),
                      conditionalPanel(
                        condition = "input.MP != 'Classification' && input.MP != 'Expression per gene'",
                        ## Which node to display 2
                        uiOutput(outputId = "nodebutton")
                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Profile scores' || input.MP == 'Genes used by CHETAH'",
                        ## Which type to display 3
                        uiOutput(outputId = "typebutton")
                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Genes used by CHETAH'",
                        ## Scale the matrix 3
                        checkboxInput(inputId = "scaling",
                                      label = p('Scale Matrix', class = 'opt'),
                                      value = FALSE),
                        ## Number of genes 3
                        numericInput(inputId = "n_genes",
                                     label = p('# of genes', class = 'opt'),
                                     value = 60, step = 1, min = 2, max = 200),
                        ## Height of the matrix 3
                        numericInput(inputId = "lettersize",
                                     label = p('Matrix label size', class = 'opt'),
                                     value = 10, step = 0.5, min = 1, max = 20)
                      ),
                      conditionalPanel(
                        condition = "input.MP == 'Genes used by CHETAH'",
                        ## Height of the matrix 3
                        checkboxInput(inputId = "largediff",
                                      label = p('Genes with max. difference in the INPUT', class = 'opt'),
                                      value = TRUE),
                        checkboxInput(inputId = "inclnodes",
                                      label = p('Add intermediate types', class = 'opt'),
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
                                           tags$div(title =
"The cells on the (t-SNE) coordinates that you provided.
Each cell is colored according to the CHETAH cell type.
As in any plot, the strictness of the classification
can be changed using 'Confidence Threshold' on the left panel.
Also, inferred intermediate types (nodes in the tree below)
can be colored with 'Color the intermediate types'.

Please note that cells labeled: 'Unknown' are cells for which
classification stopped in Node0

If more cell types are present than could be shown,
a small scroll bar is present on the right side of the plot",
                                              h3('CHETAH classification')
                                           ),
                                           plotly::plotlyOutput(outputId = 'classTsne', height = 600),
                                           downloadButton('dwn_clTsne', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Classification percentages
                                           tags$div(title =
"The percentages of cell types in you sample,
according to the current 'Confidence Threshold'",
                                               h3("Cell type percentages")
                                           ),
                                           plotOutput('typestable'),
                                           downloadButton('dwn_ttable', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Classification tree
                                           tags$div(title =
"The classification tree that CHETAH constructed
by hierarchically clustering the reference data.
The (leaf) nodes are colored with the same colors
as used in the plots above.",
                                              h3("Classification tree")
                                           ),
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
                                           tags$div(title =
"Same as on previous page, but the cells are colored
by the confidence score of the currently selected node.
If a cell is assigned to the left branch,
the scores are negative (for plotting purposes)
and the colors shades of blue.
For the right branch shades of red are used.",
                                              h3('Confidence scores per nodes')
                                           ),
                                           plotly::plotlyOutput(outputId = "confTsne", height = 600),
                                           downloadButton('dwn_confTsne', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(4,
                                           tags$div(title =
"The classification tree.
Blue cells in the t-SNE plot above are assigned to the blue branches,
the red cells are assigned to the red branches.
The branch that is shown can be changed by 'Choose Node'",
                                              h3('Current branches')
                                           ),
                                           ## Colored cl. tree
                                           plotOutput(outputId = "confTree", height = 600),
                                           downloadButton('dwn_confTree', 'Download Plot')
                                    ),
                                    column(8,
                                           tags$div(title =
"A heatmap of the Profile scores that are used
for the Confidence score calculation.
For an input cell (column), the more strongly the colors/scores differ
between reference types of the two branches,
the more confidence a assignment is.
The branch that is shown can be changed by 'Choose Node'",
                                              h3('Profile score heatmap')
                                           ),
                                           ## Heatmap
                                           plotOutput(outputId = "confHM", height = 600),
                                           tags$div(title =
"Check to sort alphabetically by cell type.
If not checked, cells are sorted by confidence score: 
From highest for the left branch, to highest for the right branch.",
                                           checkboxInput(inputId = "sortByType",
                                                         label = p('Sort by Type', class = 'opt'),
                                                         value = FALSE)),
                                           downloadButton('dwn_confHM', 'Download Plot')
                                    )
                                  )
                         ),
                         ### Tab 3
                         tabPanel("Profile scores",
                                  fluidRow(
                                    column(12,
                                           ## Profile score plots
                                           tags$div( title =
"Same as in the previous tabs,
but the cells are colored by the profile score
of the currently selected reference profile.
Change the selected type by adjusting
the 'Choose Node' and
'Choose Type' buttons",
                                              h3('Profiles scores per nodes per type')
                                           ),
                                           plotly::plotlyOutput(outputId = "profTsne", height = 600),
                                           downloadButton('dwn_profTsne', 'Download Plot')
                                           )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Profile scores in boxplots
                                           tags$div(title =
"The same values as above, but in box plots, grouped by the CHETAH cell type.",
                                              h3('In a boxplot')
                                           ),
                                           plotOutput(outputId = "profBox", height = 600),
                                           downloadButton('dwn_profBox', 'Download Plot')
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           ## Tree
                                           tags$div(title =
"The tree, with only the selected reference profile
colored (red if right branch, blue if left)
and the types with which the selected type is compared have the opposite color.
Is a Profile Score high for type X (with Y and Z in the other branch),
than this means that this cell is more likely
to be of type X than of types Y and Z.",
                                              h3('Tree')
                                           ),
                                           plotOutput(outputId = "profTree", height = 600)
                                    )
                                  )
                         ),
                         ### Tab 4
                         tabPanel("Genes used by CHETAH",
                                  fluidRow(
                                    column(12,
                                           ## Gene heatmap
                                           tags$div(title =
"The expression of the genes that were used by CHETAH
for the classification in a heatmap. ,
CHETAH selects a differnt group of genes
for each reference type, in each node.
The cell types to which the current type is compared,
can be viewed in the tree below the heatmap.
Change the selected type by adjusting
the 'Choose Node' and 'Choose Type' buttons
The matrix can be scaled per row (gene), by clicking 'Scale Matrix'.
The number of genes to show can be adjusted by '# of genes' (max = 200).
Automatically, the genes that had the highest difference in the input data will be selected.
To select the genes with the highest difference in the reference data:
check the 'Genes with max. difference in input' box
Both the plot and a list of the selected genes can be downloaded
(buttons below the plot)",
                                              h3("Selected genes per cell type")
                                           ),
                                           plotOutput(outputId = "exprHM", height = 800),
                                           downloadButton('dwn_exprHM', 'Download Plot'),
                                           downloadButton('dwn_genes', 'Download Genes'),
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
                                             tags$div(title =
"Select a gene, for which to plot 
the expression per cell, grouped
and colored per cell type in boxplots",
                                                h3("Gene expression per cell type")
                                             ),
                                             plotOutput(outputId = "geneExpr", height = 400),
                                             downloadButton('dwn_geneExpr', 'Download Plot')
                                      )
                                  )
                          ),
                          ### Tab 6
                          tabPanel("Info",
                                   fluidRow(
                                      column(2),
                                      column(10,
                                         tags$div(style = "right:100px;width:600px;text-align:justify;",
                                           tags$h3('The general method'),
                                           tags$p("CHETAH creates a classifcation tree, by hierarchically clustering the reference data.
                                           CHETAh calculates two types of scores for each input cell j in each node of the classification tree.
                                           First, for each cell j, a Profile Score is calculated for each reference cell type in this node
                                           This Profile Score is a -1:1 value that represents the similarity of j to the reference profile.
                                           Based of these Profile Scores, cell j is assigned to the right or the left branch of the node.
                                           For this assignment, one Confidence Score is calculated. This score (values of 0:2) represents the confidence of the assignment
                                           If the assignment of j has a confidence below the Confidence Threshold, the cell is not assigned to a branch, but for j, the classification will stop in that node.
                                           A classifcation to a non-leaf node are called an 'intermediate type'.
                                           Classifcations to a leaf node of the tree is called a 'final type'.
                                           0.1 is the default Confidence threshold. A threshold of 0 will classify all cells to a final type (leaf node of the tree) and 1 is a very stringent cut-off,
                                           in which case only the cells with very high confidence will be classified."),
                                           tags$br(),
                                           tags$p("Recap: for each cell j, in each node, a Profile Score is calculated for each reference profile, the cell is assigned to a branch, and for this assignment, one Confidence Score is calculated."),
                                           tags$br(),
                                           tags$p("Most plots have a help-message that can be viewed by hovering over the title of the plot.")
                                         )
                                      )
                                   )
                          )
                )
            )
        )
    )
)

# Define UI ----
ui <- fluidPage(
  tags$head(
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
                      checkboxInput(inputId = "colorsplits",
                                    label = p('Color the intermediate types', class = 'opt'),
                                    value = FALSE),
                      ## Conf threshhold for which cells to classify
                      numericInput("conf_thresh", label = p("Confidence threshold", class = 'opt'),
                                   value = 0.1, step = 0.01),

                      conditionalPanel(
                        condition = "input.MP != 'Genes used by CHETAH' && input.MP != 'Differentially expressed genes'",
                        ## t-SNE point size
                        numericInput("ptsize", label = p("Point Size", class = 'opt'),
                                     value = 1.5, step = 0.5)
                      ),
                      conditionalPanel(
                        condition = "input.MP != 'Classification' && input.MP != 'Differentially expressed genes'",
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
                                           plotOutput('typestable')
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
                         tabPanel("Differentially expressed genes",
                                  fluidRow(
                                    column(12,
                                           ## Gene heatmap
                                           h3("Differentially expressed genes per cell type"),
                                           plotOutput(outputId = "diff_exp_HM", height = 800),
                                           downloadButton('dwn_diff_expr_HM', 'Download Plot')
                                    )
                                  )
                         )
             )
           )
    )
  )
)

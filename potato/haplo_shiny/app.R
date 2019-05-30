packages <- c("ggplot2", "shiny", "igraph", "devtools", "stringr", "tidyverse", "grid", "ggplotify")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
if(!require(mappoly))
{
  install_github("mmollina/mappoly")  
  require(mappoly)
}
require(igraph)
require(shiny)
require(ggplot2)
require(stringr)
require(tidyverse)
require(gtable)
require(grid)
require(mappoly)
source("detect_meiotic_configuration_functions.R")
source("global.R")
# Example of UI with fluidPage
ui <- fluidPage(
  # Application title
  titlePanel("Haplotype reconstruction and multivalent evidence"),
  sidebarLayout(
    sidebarPanel(
      selectInput('ch', 'Chromosome', c(1:12)),
      selectInput('ind', 'Individual', ind.names),
      sliderInput("hom.prob.thresh", "Homologous probability threshold",
                  min = 0.5, max = 1,
                  value = 0.80, step = 0.05),
      sliderInput("seg.length.thresh", "Segment/gap length threshold (cM)",
                  min = 0, max = 30,
                  value = 5, step = 1),
      sliderInput("dist.thresh", "Minumum distance to differentiate c.o. (cM)",
                  min = 0, max = 2,
                  value = .5, step = 0.1),
      sliderInput("perc.info", "Minimum percentage of informative regions",
                  min = 0, max = 100,
                  value = 20, step = 10),
      checkboxInput("op", label = "Show original probability profile", value = FALSE),
      includeHTML("include2.html")
    ),
    mainPanel(h4("Homologous probability profiles"),
              fluidRow(
                splitLayout(cellWidths = c("49.5%", "49.5%","1%"), HTML("P1"), HTML("P2"),""),
                splitLayout(cellWidths = c("50%", "50%"), conditionalPanel("input.op", plotOutput('plotgraph01', height = "250px")),
                            conditionalPanel("input.op", plotOutput('plotgraph02', height = "250px"))),
                splitLayout(cellWidths = c("50%", "50%"), plotOutput("plotgraph1", height = "250px"), plotOutput("plotgraph2", height = "250px")),
                h4("Representation of the meiotic results"),
                splitLayout(cellWidths = c("50%", "50%"), plotOutput("plotgraph001", height = "200px"), plotOutput("plotgraph002", height = "200px")),
                HTML("The graph nodes represent homologous chromosomes and edges represent recombination events between them."),
                h4("Meiosis summary"),
                splitLayout(cellWidths = c("49.5%", "49.5%","1%"), verbatimTextOutput("txtout1"), verbatimTextOutput("txtout2"),""),
                splitLayout(cellWidths = c("49.5%", "49.5%","1%"), htmlOutput("text_b"), htmlOutput("text_t"),""),
                h4("Notes"),
                includeHTML("include.html")
              )
    )
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize=300*1024^2) 
  selectedData <- reactive({
    load(file = paste0("./dat/hom_prob_ch", input$ch, ".rda"))
    return(list(zb = round(zb[input$ind,,],2), zt = round(zt[input$ind,,],2), map = map))
  })
  loadMap<-reactive({
    load(paste0("./dat/dat_ind_",input$ind,".rda"))
    load(paste0("./dat/map_ch",input$ch,".rda"))
    return(list(map = map.mappoly, dat = dat.mappoly))
  })
  plot_homologous<-reactive({
    x1<-plot_recombination_points(pr.hom = selectedData()$zb, 
                              map = selectedData()$map, title.plot = "",
                              individual = input$ind, 
                              parent = "P1",
                              hom.prob.thresh = input$hom.prob.thresh, 
                              seg.length.thresh = input$seg.length.thresh, 
                              perc.info = input$perc.info,
                              map.mappoly = loadMap()$map, dat.mappoly = loadMap()$dat,
                              dist.thresh = input$dist.thresh)
    x2<-plot_recombination_points(pr.hom = selectedData()$zt, 
                                  map = selectedData()$map, 
                                  individual = input$ind, 
                                  parent = "P2",  title.plot = "",
                                  hom.prob.thresh = input$hom.prob.thresh, 
                                  seg.length.thresh = input$seg.length.thresh, 
                                  perc.info = input$perc.info,
                                  map.mappoly = loadMap()$map, dat.mappoly = loadMap()$dat,
                                  dist.thresh = input$dist.thresh)
    return(list(original.P1 = x1$plot.orig, original.P2 = x2$plot.orig,
                plot.P1 = x1$plot, summary.P1 = x1$summary, 
                plot.P2 = x2$plot, summary.P2 = x2$summary,
                meiotic.configuration.b = x1$meiotic.configuration, 
                meiotic.configuration.t = x2$meiotic.configuration, 
                meiotic.graph.b = x1$meiotic.graph, meiotic.graph.t = x2$meiotic.graph))
  })
  output$plotgraph001 <- renderPlot({ par(mai = c(0,0,0,0))
                                         plot(plot_homologous()$meiotic.graph.b$g,
                                         edge.arrow.size=.2,
                                         vertex.label.color="black",
                                         vertex.color=plot_homologous()$meiotic.graph.b$vertex.color)})
  output$plotgraph002 <- renderPlot({ par(mai = c(0,0,0,0))
                                      plot(plot_homologous()$meiotic.graph.t$g,
                                         edge.arrow.size=.2,
                                         vertex.label.color="black",
                                         vertex.color=plot_homologous()$meiotic.graph.t$vertex.color)})
  output$plotgraph01 <- renderPlot(plot_homologous()$original.P1)
  output$plotgraph02 <- renderPlot(plot_homologous()$original.P2)
  output$text_b <- renderText({ 
    if(plot_homologous()$meiotic.configuration.b == "t"){
      return(paste('<button type="button" class="btn btn-labeled btn-warning"><span class="btn-label"><i class="glyphicon glyphicon-ok"></i></span>  Evidence for quadrivalent formation.</button>'))}
    else if(plot_homologous()$meiotic.configuration.b == "h"){
      return(paste('<button type="button" class="btn btn-labeled btn-danger"><span class="btn-label"><i class="glyphicon glyphicon-ok"></i></span>  Evidence for hexavalent formation.</button>'))} 
    else{return(paste('<button type="button" class="btn btn-labeled btn-info"><span class="btn-label"><i class="glyphicon glyphicon-remove"></i></span>  No evidence of multivalent formation.</button>'))}})
  output$text_t <- renderText({ 
    if(plot_homologous()$meiotic.configuration.t == "t"){
      return(paste('<button type="button" class="btn btn-labeled btn-warning"><span class="btn-label"><i class="glyphicon glyphicon-ok"></i></span>  Evidence for quadrivalent formation.</button>'))}
    else if(plot_homologous()$meiotic.configuration.t == "h"){
      return(paste('<button type="button" class="btn btn-labeled btn-danger"><span class="btn-label"><i class="glyphicon glyphicon-ok"></i></span>  Evidence for hexavalent formation.</button>'))} 
    else{return(paste('<button type="button" class="btn btn-labeled btn-info"><span class="btn-label"><i class="glyphicon glyphicon-remove"></i></span>  No evidence of multivalent formation.</button>'))}})
  output$plotgraph1 <- renderPlot(plot_homologous()$plot.P1)
  output$plotgraph2 <- renderPlot(plot_homologous()$plot.P2)
  output$txtout1 <- renderPrint({plot_homologous()$summary.P1})
  output$txtout2 <- renderPrint({plot_homologous()$summary.P2})
}

shinyApp(ui, server)
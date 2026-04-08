#' Shiny app of ccTools
#'
#' Run this function to interactively run ccTools in a shiny app
#' @export
ccTools<-function(){
  ui<-shiny::bootstrapPage(
    shiny::titlePanel("ccTools"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput(inputId="datfile",
                         label="Select File:",
                         accept=".txt"),
        shiny::uiOutput('refGen'),
        shiny::uiOutput('naiveLev'),
        shiny::radioButtons(inputId="rmOut",
                            label="Remove Outliers?",
                            choices=c("Yes","No"),
                            selected="Yes"),
        shiny::numericInput(inputId="FDR",
                            label="FDR threshold for outliers: (ROUT method):",
                            min=0,
                            max=1,
                            value=0.05)
      ),
      shiny::mainPanel(
        shiny::fluidRow(
          shiny::column(12, shiny::plotOutput(outputId="boxplot"))
        )
      )
    )
  )

  server<-function(input, output){
    filename<-shiny::reactive(
      return(input$datfile$datapath)
    )
    getGen<-shiny::reactive({
      fn<-filename()
      x<-ccTools::loadData(datname=fn, refgenotype="return")
      return(x)
    })
    output$refGen<-shiny::renderUI({
      opt<-getGen()
      selectInput(inputId="genRef",
                  label="Reference Genotype:",
                  multiple=F,
                  choices=opt,
                  selected=opt[1])
    })
    getLev<-shiny::reactive({
      fn<-filename()
      x<-ccTools::loadData(datname=fn, naivelevel = "return")
      return(x)
    })
    output$naiveLev<-shiny::renderUI({
      opt<-getLev()
      selectInput(inputId="levNaive",
                  label="Naive Level:",
                  multiple=F,
                  choices=opt,
                  selected=opt[1])
    })
    output$boxplot<-shiny::renderPlot({
      out<-F
      if(input$rmOut=="No"){
        out<-T
      }
      fn<-filename()
      dat<-ccTools::loadData(datname=fn, refgenotype = input$genRef, naivelevel = input$levNaive)
      dat<-ccTools::findOutlier(dat, Q=input$FDR, outliers=out)
      ccBoxplots(dat)
    })
  }

  shiny::shinyApp(ui=ui, server=server)
}

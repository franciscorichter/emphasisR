rm(brts)

ui <- fluidPage(
  sidebarLayout(position = "left",
                sidebarPanel("Controls",
                             textInput('vec1', 'Enter a vector (comma delimited) with branching times', "0,1,2"),
                             actionButton("gogobutt","Go"),
                             actionButton("stopbutt","Stop"),
                             actionButton("resetbutt","Reset"),
                             numericInput("ss", "Monte-Carlo sample size:", 10),
                             #numericInput("brts", "brts:", c(0.1,0.2,3,4)),
                             numericInput("lambda", "Initial lambda:", 10)),
  
                mainPanel("Plot",
                          plotOutput("histplot"),
                          plotOutput("mu"),
                          plotOutput("K")
                          )
  ))
server <- function(input,output,session) {
  pars = c(10,0.1,5)
  rv <- reactiveValues(x=c(10,0.1,5),run=F)
  autoInvalidate <- reactiveTimer(intervalMs=500,session)
  
  observe({
    autoInvalidate()
    isolate({ if (rv$run) { 
      brts <- as.numeric(unlist(strsplit(input$vec1,",")))
      rv$x <- rbind(rv$x,mcem_step(brts,pars,maxnumspec = 15,MC_ss = input$ss))
      pars = rv$x[nrow(rv$x),]
    } })
  })
  
  observeEvent(input$gogobutt, { isolate({ rv$run=T      }) })
  observeEvent(input$stopbutt, { isolate({ rv$run=F      }) })
  observeEvent(input$resetbutt,{ isolate({ rv$x=mcem_step(brts,c(5,0.2,10),maxnumspec = 15,MC_ss = input$ss) }) })
  
  output$histplot <- renderPlot({
    htit <- sprintf("Hist of %d rnorms",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,1],type="l")
    ##hist(rv$x,col = "steelblue",main=htit,breaks=12)
  })
  output$mu <- renderPlot({
    htit <- sprintf("Hist of %d rnorms",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,2],type="l")
    ##hist(rv$x,col = "steelblue",main=htit,breaks=12)
  })
  output$K <- renderPlot({
    htit <- sprintf("K",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,3],type="l")
  })
}
shinyApp(ui, server)
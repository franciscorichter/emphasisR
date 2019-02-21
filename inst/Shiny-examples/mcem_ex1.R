ui <- fluidPage(
  sidebarLayout(position = "left",
                sidebarPanel("Controls",
                             actionButton("gogobutt","Go"),
                             actionButton("stopbutt","Stop"),
                             actionButton("resetbutt","Reset"),
                             numericInput("ss", "Monte-Carlo sample size:", 10),
                             numericInput("lambda", "Initial lambda:", 10)),
                mainPanel("Plot",
                          plotOutput("histplot"),
                          plotOutput("mu"))
  ))
server <- function(input,output,session) {
  pars = c(10,0.1,5)
  rv <- reactiveValues(x=c(10,0.1,5),run=F)
  
  autoInvalidate <- reactiveTimer(intervalMs=500,session)
  
  observe({
    autoInvalidate()
    isolate({ if (rv$run) { rv$x <- rbind(rv$x,mcem_step(brts,pars,maxnumspec = 15,MC_ss = input$ss))
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
}
shinyApp(ui, server)
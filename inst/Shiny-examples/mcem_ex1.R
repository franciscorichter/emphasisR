#brts = c(0.1,0.2,3,4)
brts_d = c(4,3.9,3.8,1)

ui <- fluidPage(
  sidebarLayout(position = "left",
                sidebarPanel("Controls",
                             textInput('vec1', 'Enter a vector (comma delimited) with branching times', "0.1,0.2,3,4"),
                             actionButton("gogobutt","Go"),
                             actionButton("stopbutt","Stop"),
                             actionButton("resetbutt","Reset"),
                             numericInput("ss", "Monte-Carlo sample size:", 10)),
                             #numericInput("brts", "brts:", c(0.1,0.2,3,4)),
                             #numericInput("la", "Initial lambda:", 10)),
  
                mainPanel("Plot",
                          plotOutput("lambda"),
                          plotOutput("mu"),
                          plotOutput("K"),
                          plotOutput("fhat")
                          )
  ))
server <- function(input,output,session) {
  init_pars = c(5,0.2,10)
  rv <- reactiveValues(x=init_pars,run=F,fhat=NULL,se=NULL,ftrue=NULL,lambda=NULL)
  autoInvalidate <- reactiveTimer(intervalMs=500,session)
  pars = init_pars
  save(pars,file="first.R")
  observe({
    autoInvalidate()
    isolate({ if (rv$run) { 
      brts <- as.numeric(unlist(strsplit(input$vec1,",")))
      load("first.R")
      rv$lambda <- c(rv$lambda,pars[1])
      mcem = mcem_step(brts,pars,maxnumspec = 15,MC_ss = input$ss)
      ftrue = exp(DDD::dd_loglik(pars1 = pars, pars2 = c(250,1,0,1,0,1),brts = brts_d,missnumspec = 0))
      pars = mcem$pars
      save(pars,file="first.R")
      fhat = mcem$fhat
      se = mcem$se
      rv$x <- rbind(rv$x,pars)
      PARS = rv$x
      save(PARS,file="dinamical.RData")
      rv$fhat = c(rv$fhat,fhat)
      rv$se = c(rv$se,se)
      rv$ftrue = c(rv$ftrue,ftrue)
      #pars = rv$x[nrow(rv$x),]
    } })
  })
  
  observeEvent(input$gogobutt, { isolate({ rv$run=T      }) })
  observeEvent(input$stopbutt, { isolate({ rv$run=F      }) })
  observeEvent(input$resetbutt,{ isolate({ rv$x=mcem_step(brts,c(50,10,100),maxnumspec = 15,MC_ss = input$ss) }) })
  
  output$lambda <- renderPlot({
    htit <- sprintf("Hist of %d rnorms",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,1],type="l")
    points(1:length(rv$lambda),rv$lambda)
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
  output$fhat <- renderPlot({
    htit <- sprintf("fhat",length(rv$fhat))
    plot(1:length(rv$ftrue),rv$ftrue,ylim=c(0,max(rv$fhat,rv$ftrue,rv$fhat+1.96*rv$se)))
    lines(1:length(rv$fhat),rv$fhat,col="blue")
    lines(1:length(rv$fhat),rv$fhat+1.96*rv$se,col="red")
    lines(1:length(rv$fhat),rv$fhat-1.96*rv$se,col="red")
    
  })
}
shinyApp(ui, server)
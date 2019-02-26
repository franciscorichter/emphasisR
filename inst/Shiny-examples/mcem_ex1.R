brts_d = brts_d=c(4.9999999998,4.806886544,4.70731246478,4.50735197578,4.37856240588,4.29594855558,4.19207515688,4.18261061218,4.11238451758,4.09640902445,3.81723693538,3.71143733895,3.48845298905,3.25729503338,3.11613886835,2.64829864145,2.63531839038,2.37990087748,1.82721570435,0.83704715535,0.64242044758,0.56121103655,0.356333544350001,0.346462849050001)

#pars = DDD:::dd_ML(brts = brts,soc=1)
#pars = c(3.021189,0.164449,24.579503 )
#wt = -diff(c(brts_d,0))
#brts = cumsum(wt)

#brts = c(0.1,0.2,3,4)
#brts_d = c(4,3.9,3.8,1)

ui <- fluidPage(
  sidebarLayout(position = "left",
                sidebarPanel("Controls",
                             textInput('vec1', 'Enter a vector (comma delimited) with branching times', "0.1,0.2,3,4"),
                             numericInput("ss", "Monte-Carlo sample size:", 10),
                             #textInput('vec2', 'Enter a vector (comma delimited) with parameters', "5,0.2,10"),
                             actionButton("gogobutt","Go"),
                             actionButton("stopbutt","Stop"),
                             actionButton("resetbutt","Reset"),
                             textOutput("txtOutput")
                             ),
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
  init_pars = c(4,0.3,40)
  rv <- reactiveValues(x=init_pars,run=F,fhat=NULL,se=NULL,ftrue=NULL,lambda=NULL,LastTime=NULL)
  autoInvalidate <- reactiveTimer(intervalMs=500,session)
  pars = init_pars
  save(pars,file="first.R")
  observe({
    autoInvalidate()
    isolate({ if (rv$run) { 
      brts <- as.numeric(unlist(strsplit(input$vec1,",")))
      load("first.R")
      rv$lambda <- c(rv$lambda,pars[1])
      time = proc.time()
      mcem = mcem_step(brts,pars,maxnumspec = 35,MC_ss = input$ss)
      rv$LastTime <- get.time(time)
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
  observeEvent(input$resetbutt,{ isolate({ rv$x=mcem_step(brts,c(50,10,100),maxnumspec = 35,MC_ss = input$ss) }) })
  output$txtOutput = renderText({
    paste0("Last iteration took: ", rv$LastTime, "\n Last parameter: ", rv$lambda[length(rv$lambda)] )
    #paste0("Last parameters:", rv$lambda[length(rv$lambda)])
  })
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
    plot(1:length(rv$ftrue),rv$ftrue,ylim=c(0,max(rv$fhat,rv$ftrue)))
    lines(1:length(rv$fhat),rv$fhat,col="blue")
    lines(1:length(rv$fhat),rv$fhat+1.96*rv$se,col="red")
    lines(1:length(rv$fhat),rv$fhat-1.96*rv$se,col="red")
    
  })
}
shinyApp(ui, server)
#brts_d = brts_d=c(4.9999999998,4.806886544,4.70731246478,4.50735197578,4.37856240588,4.29594855558,4.19207515688,4.18261061218,4.11238451758,4.09640902445,3.81723693538,3.71143733895,3.48845298905,3.25729503338,3.11613886835,2.64829864145,2.63531839038,2.37990087748,1.82721570435,0.83704715535,0.64242044758,0.56121103655,0.356333544350001,0.346462849050001)
#0.1931135, 0.2926875, 0.4926480, 0.6214376, 0.7040514, 0.8079248, 0.8173894, 0.8876155, 0.9035910, 1.1827631, 1.2885627, 1.5115470, 1.7427050, 1.8838611, 2.3517014, 2.3646816, 2.6200991, 3.1727843, 4.1629528, 4.3575796, 4.4387890, 4.6436665, 4.6535372, 5.0000000


ui <- fluidPage(
  sidebarLayout(position = "left",
                sidebarPanel("Controls",
                             #textInput('vec1', 'Enter a vector (comma delimited) with branching times', "0.1,0.2,3,4"),
                             selectInput("brts", "Choose Phylo/Branching times:",
                                         list("Dendroica" = "0.1931135, 0.2926875, 0.4926480, 0.6214376, 0.7040514, 0.8079248, 0.8173894, 0.8876155, 0.9035910, 1.1827631, 1.2885627, 1.5115470, 1.7427050, 1.8838611, 2.3517014, 2.3646816, 2.6200991, 3.1727843, 4.1629528, 4.3575796, 4.4387890, 4.6436665, 4.6535372, 5",
                                              "Simple" = "0.1,0.2,3,4")
                             ),
                             numericInput("ss", "Monte-Carlo sample size:", 100),
                             numericInput("Bt", "Number of best trees to take:", 10),
                             #textInput('vec2', 'Enter a vector (comma delimited) with parameters', "5,0.2,10"),
                             actionButton("gogobutt","Go"),
                             actionButton("stopbutt","Stop"),
                             actionButton("resetbutt","Reset"),
                             textOutput("txtOutput1"),
                             textOutput("txtOutput2"),
                             textOutput("txtOutput3")
                             ),
                             #numericInput("brts", "brts:", c(0.1,0.2,3,4)),
                             #numericInput("la", "Initial lambda:", 10)),
          
                mainPanel("Parameters",
                          plotOutput("lambda"),
                          plotOutput("mu"),
                          plotOutput("K"),
                          "Diagnostics",
                          plotOutput("fhat"),
                          plotOutput("rellik")
                          ),
                
                
  ))
server <- function(input,output,session) {
  init_pars = c(4,0.3,40)
  rv <- reactiveValues(x=init_pars,run=F,fhat=NULL,se=NULL,ftrue=NULL,lambda=NULL,LastTime=NULL,rellik=NULL,ll.prop=NULL)
  autoInvalidate <- reactiveTimer(intervalMs=500,session)
  pars = init_pars
  save(pars,file="first.R")
  observe({
    autoInvalidate()
    isolate({ if (rv$run) { 
      brts <- as.numeric(unlist(strsplit(input$brts,",")))
      load("first.R")
      rv$lambda <- c(rv$lambda,pars[1])
      time = proc.time()
      mcem = mcem_step(brts,pars,maxnumspec = 35,MC_ss = input$ss,selectBestTrees = TRUE,bestTrees = input$Bt)
      #h1 = try(diag(solve(mcem$hessian))/m)
      rv$ll.prop = mcem$loglik.proportion
      rv$rellik = c(rv$rellik,rel.llik(S1=mcem$st$trees,p0=pars,p1=mcem$pars))
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
  observeEvent(input$resetbutt,{ isolate({ rv$x=mcem_step(as.numeric(unlist(strsplit(input$vec1,","))),c(50,10,100),maxnumspec = 35,MC_ss = input$ss) }) })
  output$txtOutput1 = renderText({
    paste0("Last iteration took: ", rv$LastTime)
  })
  output$txtOutput2 = renderText({
    paste0("Last iteration took: ", "la: ", rv$x[nrow(rv$x),1], " mu: ", rv$x[nrow(rv$x),2], " K: ", rv$x[nrow(rv$x),3])
  })
  output$txtOutput3 = renderText({
    paste0("Proportion of likelihood: ", rv$ll.prop )
  })
  output$lambda <- renderPlot({
    htit <- sprintf("Hist of %d rnorms",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,1],type="l")
    #points(1:length(rv$lambda),rv$lambda)
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
  output$rellik <- renderPlot({
    htit <- sprintf("rellik",length(rv$fhat))
    plot(1:length(rv$rellik),rv$rellik,type="l")
    
  })
}
shinyApp(ui, server)
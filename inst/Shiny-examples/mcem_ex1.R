#brts_d = brts_d=c(4.9999999998,4.806886544,4.70731246478,4.50735197578,4.37856240588,4.29594855558,4.19207515688,4.18261061218,4.11238451758,4.09640902445,3.81723693538,3.71143733895,3.48845298905,3.25729503338,3.11613886835,2.64829864145,2.63531839038,2.37990087748,1.82721570435,0.83704715535,0.64242044758,0.56121103655,0.356333544350001,0.346462849050001)
#0.1931135, 0.2926875, 0.4926480, 0.6214376, 0.7040514, 0.8079248, 0.8173894, 0.8876155, 0.9035910, 1.1827631, 1.2885627, 1.5115470, 1.7427050, 1.8838611, 2.3517014, 2.3646816, 2.6200991, 3.1727843, 4.1629528, 4.3575796, 4.4387890, 4.6436665, 4.6535372, 5.0000000
#plethodon
brts_d = c(11.3,9.55365380008198,9.26434040327225,8.83592350352767,8.3434446982257,8.17781491491496,7.95978190384214,6.61207494082374,6.5679856688767,6.21838471981418,5.59809615547134,5.37012669852355,4.7638222125791,4.10749650972075,4.02367324807484,3.65931960175062,3.32916292401100,3.23132222435799,3.18206288699248,2.8572235287017,2.58222342582278,2.43078192215161,1.87377417677032,1.79734091086791,1.77566693721338,1.52675067868777,1.11116172787207,0.800771123741394,0.498973146096477)

ui <- fluidPage(
  sidebarLayout(position = "left",
                sidebarPanel("Controls",
                             #textInput('vec1', 'Enter a vector (comma delimited) with branching times', "0.1,0.2,3,4"),
                             selectInput("brts", "Choose Phylo/Branching times:",
                                         list("Dendroica" = "0.1931135, 0.2926875, 0.4926480, 0.6214376, 0.7040514, 0.8079248, 0.8173894, 0.8876155, 0.9035910, 1.1827631, 1.2885627, 1.5115470, 1.7427050, 1.8838611, 2.3517014, 2.3646816, 2.6200991, 3.1727843, 4.1629528, 4.3575796, 4.4387890, 4.6436665, 4.6535372, 5",
                                              "Anolis" = "103.31057277,97.96837055,94.923127866,94.866216796,90.810203432,90.44410561,90.080885176,86.938065219,83.192481566,79.903508082,78.144291981,75.916079896,75.270062039,74.19732113,72.825735377,72.5234868711,68.360444962,64.1335159681,63.557121926,63.523319671,63.398403586,60.3209181541,59.490443993,57.576137962,56.933279789,56.6964480574,56.361043545,55.506578472,54.495983232,54.199692129,54.051109589,53.692672898,52.897212802,52.23186871,52.164805405,51.94779018,50.553819579,48.373129976,47.4174457904,47.189946167,46.7942740811,44.287638517,43.296282982,43.242616701,42.272773859,42.1266648041,41.905974158,41.329061036,41.1257958974,39.697767108,39.677765636,39.397083778,37.911582502,37.397487349,34.605557178,32.824929892,32.228763421,31.561908554,30.308206955,30.281651159,30.1639183904,29.8042173411,29.773786118,29.6447104204,29.541373926,29.5407793691,28.623740578,28.506256108,27.105138539,26.439039467,26.378922306,26.285718299,26.233564599,24.159222149,22.974342026,21.370573573,21.247374251,20.091077714,19.6384466391,19.483057152,18.985702599,16.272845218,15.9237660074,15.8035460904,15.7840819411,15.2801426021,15.08030346650,13.1506065141,12.572643169,11.2235612480000,10.5640701490000,10.30031894,9.651984572,9.577378633,9.50317724299997,8.43806557499997,6.837926447,5.73566929399999,5.38425275100001,4.4685701345,4.33572126899998,0.828930356499995,0.552543471999996",
                                              "Cetacea" = "35.857845,33.799004,32.390661,31.621529,28.000001,26.063017,26.000001,24.698215,22.044392,19.195664,18.226421,18.023412,17.939427,17.890656,16.066686,15.669702,15.099325,14.540283,14.061555,13.042869,12.847396,11.382959,11.079292,11.028304,10.70277,10.472235,9.438013,8.925942,8.81602,8.803019,8.716039,8.252102,8.20904999999999,8.143058,8.100266,7.677423,7.514394,7.176679,6.975185,6.37563,6.28945,6.04702999999999,5.897506,5.796585,5.616381,5.49323999999999,5.466433,5.27807199999999,5.26532599999999,5.263495,5.096383,4.985876,4.947171,4.927157,4.732447,4.57089299999999,4.45271899999999,4.35571699999999,4.32202299999999,4.170967,4.166225,4.045704,3.791853,3.70627600000000,3.62611499999999,3.44535999999999,3.29116399999999,3.21256099999999,3.07916999999999,3.04867399999999,2.919779,2.83297999999999,2.19441299999999,2.09621900000000,1.93481199999999,1.82123899999999,1.622022,1.570433,1.50673699999999,1.47078099999999,1.36135800000000,1.26811400000000,1.01116300000000,0.924861999999997,0.347030000000004,0.283069999999995",
                                              "Heliconius" = "16.761439,15.160156,14.40605,13.815308,13.486476,13.164424,12.373104,10.840648,10.142988,9.911296,9.273213,9.264573,9.142266,8.536825,8.441098,8.17086,7.92524,7.478269,7.255542,6.851304,5.335066,5.335061,5.152996,4.643518,4.506785,4.446959,3.780976,3.768737,3.488772,3.398945,2.433015,2.048552,1.930075,1.602332,1.302335,1.03376100000000,0.884346",
                                              "Plethodon" = "11.3,9.55365380008198,9.26434040327225,8.83592350352767,8.3434446982257,8.17781491491496,7.95978190384214,6.61207494082374,6.5679856688767,6.21838471981418,5.59809615547134,5.37012669852355,4.7638222125791,4.10749650972075,4.02367324807484,3.65931960175062,3.32916292401100,3.23132222435799,3.18206288699248,2.8572235287017,2.58222342582278,2.43078192215161,1.87377417677032,1.79734091086791,1.77566693721338,1.52675067868777,1.11116172787207,0.800771123741394,0.498973146096477",
                                              
                                              "Simple" = "0.1,0.2,3,4")
                             ),
                             numericInput("ss", "Monte-Carlo sample size:", 100),
                             numericInput("Bt", "Number of best trees to take:", 10),
                             numericInput("par1", "Initial lambda:", 1),
                             numericInput("par2", "Initial mu:", 0.1),
                             numericInput("par3", "Initial K:", 40),
                             #textInput('vec2', 'Enter a vector (comma delimited) with parameters', "5,0.2,10"),
                             actionButton("gogobutt","Go"),
                             actionButton("stopbutt","Stop"),
                             actionButton("resetbutt","Reset"),
                             textOutput("txtOutput1"),
                             checkboxInput("ddd", "Compare with DDD", FALSE),
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
  #init_pars = c(input$par1,input$par2,input$par3)
  init_pars = c(1,0.1,40)
  rv <- reactiveValues(x=init_pars,run=F,fhat=NULL,se=NULL,ftrue=NULL,lambda=NULL,LastTime=NULL,rellik=NULL,ll.prop=NULL)
  autoInvalidate <- reactiveTimer(intervalMs=500,session)
  pars = init_pars
  save(pars,file="first.R")
  observe({
    autoInvalidate()
    isolate({ if (rv$run) { 
      brts <- as.numeric(unlist(strsplit(input$brts,",")))
      if(max(brts)==brts[1]){
        wt = -diff(c(brts,0))
        brts = cumsum(c(wt))
      }
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
    paste0( "la: ", rv$x[nrow(rv$x),1], " mu: ", rv$x[nrow(rv$x),2], " K: ", rv$x[nrow(rv$x),3])
  })
  output$txtOutput3 = renderText({
    paste0("Proportion of likelihood: ", rv$ll.prop )
  })
  output$lambda <- renderPlot({
    htit <- sprintf("Hist of %d rnorms",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,1],type="l")
  #  abline(b = 0,a = 0.523708)
    #points(1:length(rv$lambda),rv$lambda)
    ##hist(rv$x,col = "steelblue",main=htit,breaks=12)
  })
  output$mu <- renderPlot({
    htit <- sprintf("Hist of %d rnorms",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,2],type="l")
  #  abline(b = 0,a = 0.025529)
    ##hist(rv$x,col = "steelblue",main=htit,breaks=12)
  })
  output$K <- renderPlot({
    htit <- sprintf("K",length(rv$x))
    plot(1:nrow(rv$x),rv$x[,3],type="l")
  #  abline(b=0,a=31.723702)
  })
  output$fhat <- renderPlot({
    htit <- sprintf("fhat",length(rv$fhat))
    plot(1:length(rv$fhat),rv$fhat,col="blue",type="l")
    if(input$ddd) points(1:length(rv$ftrue),rv$ftrue)
    lines(1:length(rv$fhat),rv$fhat+1.96*rv$se,col="red")
    lines(1:length(rv$fhat),rv$fhat-1.96*rv$se,col="red")
    
  })
  output$rellik <- renderPlot({
    htit <- sprintf("rellik",length(rv$fhat))
    plot(1:length(rv$rellik),rv$rellik,type="l")
    
  })
}
shinyApp(ui, server)
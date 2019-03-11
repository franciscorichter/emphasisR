
    

if (file.exists("first.R")) file.remove("first.R")

#require("RCurl")
#download.file("https://github.com/franciscorichter/emphasis/blob/master/inst/Shiny-examples/BirdTree.tre",destfile = "birds.tre")
phy = read.nexus(file="BirdTree.tre")
brts_birds = sort(branching.times(phy))
n_cores = detectCores()
rv = NULL
theme_emphasis =  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ui <- fluidPage(
  tags$head(tags$style(
    HTML('
         #sidebar {
         background-color:  #D8FCFC;
         }
         
         body, label, input, button, select { 
         font-family: "Arial";
         }')
  )),
  sidebarLayout(position = "left",
               
                sidebarPanel(id="sidebar",
                  img(src='logo.png', align = "center",height = 170),
                            h3("Controls"),
                             actionButton("gogobutt","Go"),
                             actionButton("stopbutt","Stop"),
                             #     actionButton("resetbutt","Reset"),
                             
                             h3("Data"),  
                             selectInput("brts", "Choose Phylo/Branching times:",
                                         list("Other" = 0,
                                              "Dendroica" = "5,4.806886544,4.70731246478,4.50735197578,4.37856240588,4.29594855558,4.19207515688,4.18261061218,4.11238451758,4.09640902445,3.81723693538,3.71143733895,3.48845298905,3.25729503338,3.11613886835,2.64829864145,2.63531839038,2.37990087748,1.82721570435,0.83704715535,0.64242044758,0.56121103655,0.356333544350001,0.346462849050001",
                                              "Anolis" = "103.31057277,97.96837055,94.923127866,94.866216796,90.810203432,90.44410561,90.080885176,86.938065219,83.192481566,79.903508082,78.144291981,75.916079896,75.270062039,74.19732113,72.825735377,72.5234868711,68.360444962,64.1335159681,63.557121926,63.523319671,63.398403586,60.3209181541,59.490443993,57.576137962,56.933279789,56.6964480574,56.361043545,55.506578472,54.495983232,54.199692129,54.051109589,53.692672898,52.897212802,52.23186871,52.164805405,51.94779018,50.553819579,48.373129976,47.4174457904,47.189946167,46.7942740811,44.287638517,43.296282982,43.242616701,42.272773859,42.1266648041,41.905974158,41.329061036,41.1257958974,39.697767108,39.677765636,39.397083778,37.911582502,37.397487349,34.605557178,32.824929892,32.228763421,31.561908554,30.308206955,30.281651159,30.1639183904,29.8042173411,29.773786118,29.6447104204,29.541373926,29.5407793691,28.623740578,28.506256108,27.105138539,26.439039467,26.378922306,26.285718299,26.233564599,24.159222149,22.974342026,21.370573573,21.247374251,20.091077714,19.6384466391,19.483057152,18.985702599,16.272845218,15.9237660074,15.8035460904,15.7840819411,15.2801426021,15.08030346650,13.1506065141,12.572643169,11.2235612480000,10.5640701490000,10.30031894,9.651984572,9.577378633,9.50317724299997,8.43806557499997,6.837926447,5.73566929399999,5.38425275100001,4.4685701345,4.33572126899998,0.828930356499995,0.552543471999996",
                                              "Cetacea" = "35.857845,33.799004,32.390661,31.621529,28.000001,26.063017,26.000001,24.698215,22.044392,19.195664,18.226421,18.023412,17.939427,17.890656,16.066686,15.669702,15.099325,14.540283,14.061555,13.042869,12.847396,11.382959,11.079292,11.028304,10.70277,10.472235,9.438013,8.925942,8.81602,8.803019,8.716039,8.252102,8.20904999999999,8.143058,8.100266,7.677423,7.514394,7.176679,6.975185,6.37563,6.28945,6.04702999999999,5.897506,5.796585,5.616381,5.49323999999999,5.466433,5.27807199999999,5.26532599999999,5.263495,5.096383,4.985876,4.947171,4.927157,4.732447,4.57089299999999,4.45271899999999,4.35571699999999,4.32202299999999,4.170967,4.166225,4.045704,3.791853,3.70627600000000,3.62611499999999,3.44535999999999,3.29116399999999,3.21256099999999,3.07916999999999,3.04867399999999,2.919779,2.83297999999999,2.19441299999999,2.09621900000000,1.93481199999999,1.82123899999999,1.622022,1.570433,1.50673699999999,1.47078099999999,1.36135800000000,1.26811400000000,1.01116300000000,0.924861999999997,0.347030000000004,0.283069999999995",
                                              "Heliconius" = "16.761439,15.160156,14.40605,13.815308,13.486476,13.164424,12.373104,10.840648,10.142988,9.911296,9.273213,9.264573,9.142266,8.536825,8.441098,8.17086,7.92524,7.478269,7.255542,6.851304,5.335066,5.335061,5.152996,4.643518,4.506785,4.446959,3.780976,3.768737,3.488772,3.398945,2.433015,2.048552,1.930075,1.602332,1.302335,1.03376100000000,0.884346",
                                              "Plethodon" = "11.3,9.55365380008198,9.26434040327225,8.83592350352767,8.3434446982257,8.17781491491496,7.95978190384214,6.61207494082374,6.5679856688767,6.21838471981418,5.59809615547134,5.37012669852355,4.7638222125791,4.10749650972075,4.02367324807484,3.65931960175062,3.32916292401100,3.23132222435799,3.18206288699248,2.8572235287017,2.58222342582278,2.43078192215161,1.87377417677032,1.79734091086791,1.77566693721338,1.52675067868777,1.11116172787207,0.800771123741394,0.498973146096477",
                                              "Foraminifera" = "64.95,64.9,61.2,44.15,37.2,30.2,23.8,22.5,21,18.8,18.3,17.3,17,14.8,10.8,10.2,10,8.2,8,7.2,5.2,4.5,3.8,3.5,3.4,3.2,2,2,0.8,0.3,0.3",
                                              "Birds" = paste(as.character(brts_birds), collapse=", "))
                                        
                             ),
                             textInput('vec1', 'Or enter a vector (comma delimited) with branching times (Selecting Other)', "4,3.9,3.8,1"),
                             h3("Settings"),
                             numericInput("ss", "Monte-Carlo sample size:", 1000),
                             numericInput("Bt", "Number of best trees to take:", 20),
                             numericInput("maxspec", "Maximum number of missing species:", 40),
                             
                             h3("Initial parameters"),
                             numericInput("par1", "Initial lambda:", 1),
                             numericInput("par2", "Initial mu:", 0.1),
                             numericInput("par3", "Initial K:", 40),
                             
                             h3("Options"),
                             checkboxInput("ddd", "Compare with DDD", FALSE),
                             conditionalPanel(length(rv$fhat)>10, checkboxInput("CI", "Check CI (after it 10)", FALSE)),
                             checkboxInput("log", "show estimated lkelihood on log scale", FALSE),
                             checkboxInput("log_w", "show weights on log scale", TRUE),
                             numericInput("charts", "See charts from iteration:", 1),
                             numericInput("cores",paste("Your computer holds",n_cores,"cores, how many of them you want to use?"),2)#,
                             
                           
                             
                ),
                
                mainPanel(
                  
                  tabsetPanel(type = "tabs",
                              tabPanel("Analysis",
                             
                                fluidRow(
                                  column(5,
                                      h3("Parameters"),
                                       plotOutput("lambda"),
                                       plotOutput("mu"),
                                       plotOutput("K")),
                                  
                                  column(5,
                                        h3("Diagnostics"),
                                        plotOutput("fhat"),
                                        plotOutput("diff_fhat_hist"),
                                        plotOutput("rellik_hist")
                                       ),
                                  column(5,
                                         h3("Weights"),
                                         plotOutput("hist_w")  #weight_vs_dimension
                                  ),
                                  column(5,
                                         h3("Information"),
                                         textOutput("txtOutput2"),
                                         textOutput("txtOutput3")
                                  )
                                )),
                              tabPanel("Help",
                                      h1("Welcome to Emphasis!"),
                                      
                                      "Emphasis diversity dependence is a branch of the Emphasis framework for a diversity dependence model.\n",
                                      "This GUI version of the frameworks allows interactive and user friendly usage for Emphasis analysis.\n",
                                      "If you are familiar with the emphasis approach for diversification you just need a short description of the functionalities of this app and we are good to go.\n",
                                      "As a simple example we suggest to type the tree 4,3.9,3.8,1. That tree has also a pre-loaded DDD solution to compare.",
                                      
                                      h3("Controls"),
                                      "To initialize and stop the MCEM routine",
                                      h3("Data"),
                                      "For diversity dependence version data consists on a vector with branching times starting from stem age in decreasing order.",
                                      'It is also possible to use one of the 7 data examples, comparison with DDD values is possible in all of them except the birds\' tree.',
                                      h3("Settings"),
                                      h4("Monte Carlo Sampling size"),
                                      "This is the Monte Carlo sampling size on the E step. You can check if it is large enough looking at the variation on the plots as well as the weights plot. It is possible to increase it at any moment of the routine, the change will apply after the current iteration next iteration.",
                                      h4("Number of best trees to take"),
                                      "Numbers of trees to take on the M-step. Check with the weights plot and the percentage of likelihood information if the current is appropriate or not.",
                                      h4("Maximum number of species"),
                                      "The current data augmentation algorithm requires a limit on the number of missing species, you can check if that limit makes sense by checking the weights plot. If the weights are spread all over the plot you probably need to increase the maximum number of species allowed to sample in the whole space. If the last meaningful weight is far to the left to your limit then you probably want to decrease this quantity for the sake of efficiency.",
                                      h3("Initial parameters"),
                                      "Initial parameters for diversity dependence model.",
                                      h3("Options"),
                                      "Several options are available for analysis",
                                      h4("Compare with DDD"),
                                      "For the given examples (except birds\' tree) a solution is given by DDD package, useful for comparison.",
                                      h4("Check confidence intervals"),
                                      "Confidence intervals for the parameters on the MCEM process can be calculated. That is considering a combination of the MC error and the EM variance. For reliable calculations of the second one we suggest to click this button after iteration 10.",
                                      h4("log of estimated likelihood"),
                                      "To observe the estimated loglikelihood",
                                      h4("log of weights"),
                                      "To observe the weights on logarithm scale",
                                      h4("Iteration to show on plot"),
                                      "To observe the plots on more detail from last iterations",
                                      h4("Number of processors"),
                                      "Number of cores to be used on parallel. We suggest 2 for small samples and increase it as the sample size increases. Ready to go to the analysis tab?"
                                      
                                      
                                       ),
                              tabPanel("About",
                                       h2(""),
                                       #img(src='logo.png', align = "left",height = 180),
                                       helpText("Emphasis is a joint work by Francisco Richter, Rampal Etienne and Ernst Wit. For questions or comments please write to f.richter@rug.nl.\n", align = "left",height = 70),
                                       img(src='license.png', align = "right",height = 60)
                              )
                  )
                )
              )
  )

            
  
              
              
server <- function(input,output,session) {
  rv <- reactiveValues(x=c(NULL,NULL,NULL),run=F,fhat=NULL,diff_fhat=NULL,se=NULL,ftrue=NULL,LastTime=NULL,rellik=NULL,ll.prop=NULL,mle_dd=c(NULL,NULL,NULL),H=c(NULL,NULL,NULL),sdl=NULL,sdm=NULL,sdk=NULL,dim=NULL,weights=NULL,logweights=NULL)
  autoInvalidate <- reactiveTimer(intervalMs=500,session)
  observe({
    MCEM = data.frame(it=0,par1=input$par1,par2=input$par2,par3=input$par3,fhat=NaN,diff.fhat=NaN,fhat.se=NaN,ftrue=NaN,times=NaN,mc.samplesize=input$ss,rel.lik=NaN,effective.size=NaN)
    init_pars = c(input$par1,input$par2,input$par3)
    pars = c(init_pars[1],init_pars[2],init_pars[3])
    autoInvalidate()
    isolate({ if (rv$run) { 
      if(input$brts==0){
        brts = brts_d <- as.numeric(unlist(strsplit(input$vec1,",")))
      }else{ 
        brts = brts_d <- as.numeric(unlist(strsplit(input$brts,",")))
      }
      if(brts[1] == 11.3) rv$mle_dd = c(0.523708,0.025529,31.723702)
      if(brts[1] == 5) rv$mle_dd = c(3.659926,0.182003,23.507341)
      if(brts[1] == 103.31057277) rv$mle_dd = c(0.068739,0.005933,110.066841)
      if(brts[1] == 35.857845) rv$mle_dd = c(0.135271,0.000160,234.978643)
      if(brts[1] == 16.761439) rv$mle_dd = c(0.457300,0.048939,37.782661)
      if(brts[1] == 4) rv$mle_dd = c(16.132848,0.102192,3.974553,-1.020150 )
      if(brts[1] == 64.95) rv$mle_dd = c(1.144462,0.115144,30.423973)
      if(max(brts)==brts[1]){
        wt = -diff(c(brts,0))
        brts = cumsum(c(wt))
      }
      if(file.exists("first.R")) load("first.R")
      time = proc.time()
      rv$x <- rbind(rv$x,pars)
      mcem = mcem_step(brts,pars,maxnumspec = input$maxspec,MC_ss = input$ss,selectBestTrees = TRUE,bestTrees = input$Bt,no_cores = input$cores)
      rv$H =  rbind(rv$H,mcem$h1/input$ss) # Hessian of the current parameters
      rv$ll.prop = mcem$loglik.proportion # proportion (on weights) of the likelihood considered for optimization
      rv$rellik = c(rv$rellik,rel.llik(S1=mcem$st$trees,p0=pars,p1=mcem$pars)) # relative lkelihood
      rv$weights = mcem$st$weights
      rv$logweights = mcem$st$logweights
      rv$dim = sapply(mcem$st$trees,FUN = function(list) length(list$to))
      rv$LastTime <- get.time(time) 
      if(length(brts_d)<800) rv$ftrue = c(rv$ftrue,exp(DDD::dd_loglik(pars1 = pars, pars2 = c(250,1,0,1,0,1),brts = brts_d,missnumspec = 0)))
      pars = mcem$pars
      save(pars,file="first.R")
      fhat = mcem$fhat
      se = mcem$se
      PARS = rv$x
     # save(PARS,file="dinamical.RData")
      rv$fhat = c(rv$fhat,fhat)
      rv$diff_fhat = diff(rv$fhat)
      rv$se = c(rv$se,se)
      rv$mcem_it = data.frame(it=1:length(rv$x[,1]),lambda=rv$x[,1],mu=rv$x[,2],K=rv$x[,3])
      if(length(rv$fhat)>9){
       
        gamLambda = gam(lambda ~ s(it), data=rv$mcem_it)
        gamMu = gam(mu ~ s(it), data=rv$mcem_it)
        gamK = gam(K ~ s(it), data=rv$mcem_it)
        rv$sdl = sqrt(-rv$H[,1]+gamLambda$sig2)#/input$ss)
        rv$sdm = sqrt(-rv$H[,2]+gamMu$sig2)#/input$ss)
        rv$sdk = sqrt(-rv$H[,3]+gamK$sig2)#/input$ss)
      }else{
        #i think i dont need this 
        rv$sdl = c(rv$sdl,0)
        rv$sdm = c(rv$sdm,0)
        rv$sdk = c(rv$sdk,0)
      }
      
    } })
  })
  
  observeEvent(input$gogobutt, { isolate({ rv$run=T      }) })
  observeEvent(input$stopbutt, { isolate({ rv$run=F      }) })
 # observeEvent(input$resetbutt,{ isolate({ rv$x=mcem_step(as.numeric(unlist(strsplit(input$vec1,","))),c(50,10,100),maxnumspec = input$maxspec,MC_ss = input$ss) }) })
  output$txtOutput1 = renderText({
   "Welcome to Emphasis"
  })
  output$txtOutput2 = renderText({
     if(is.null(rv$LastTime)){
    "Initializing process"
      }else{
      timescale = " sec"
      time_s = rv$LastTime
      if(rv$LastTime>60 & rv$LastTime < 3600){
        time_s = rv$LastTime/60
        timescale = " min"
      }
      if(rv$LastTime > 3600){
        time_s = rv$LastTime/3600
        timescale = " hour"
      } 
      paste0("Last iteration took: ", time_s, timescale)
      }
  })
  
  output$txtOutput3 = renderText({
    paste0("Proportion of likelihood: ", rv$ll.prop )
  })
  
  output$lambda <- renderPlot({
    if(length(rv$x)>5){ 
      if(length(rv$fhat)>9) rv$mcem_it$sdl = rv$sdl
        gl = ggplot(rv$mcem_it) + geom_line(aes(it,lambda)) + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$lambda[input$charts:length(rv$mcem_it$lambda)])),subtitle =   paste("number of last iterations to consider: ", as.character(rv$mle_dd))) 
        if(input$ddd) gl = gl + geom_hline(yintercept = rv$mle_dd[1])
        if(input$CI) gl = gl + geom_errorbar(aes(x=it, y=lambda, ymin = lambda-1.96*sdl, ymax = lambda + 1.96*sdl), colour='darkgreen') 
        gl + theme_emphasis+ geom_hline(yintercept = mean(rv$mcem_it$lambda[input$charts:length(rv$mcem_it$lambda)]), colour="blue")
        
    }
  })
  
  output$mu <- renderPlot({
    if(length(rv$x)>5){ 
      if(length(rv$fhat)>9) rv$mcem_it$sdm = rv$sdm
      gl = ggplot(rv$mcem_it) + geom_line(aes(it,mu)) + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$mu[input$charts:length(rv$mcem_it$mu)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$mu)-input$charts)) 
      if(input$ddd) gl = gl + geom_hline(yintercept = rv$mle_dd[2])
      if(input$CI) gl = gl + geom_errorbar(aes(x=it, y=mu, ymin = mu-1.96*sdl, ymax = mu + 1.96*sdl), colour='darkgreen') 
      gl + theme_emphasis + geom_hline(yintercept = mean(rv$mcem_it$mu[input$charts:length(rv$mcem_it$mu)]), colour="blue")
      
    }
  })
  
  output$K <- renderPlot({
    if(length(rv$x)>5){ 
      if(length(rv$fhat)>9) rv$mcem_it$sdk = rv$sdk
      gl = ggplot(rv$mcem_it) + geom_line(aes(it,K)) + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      if(input$ddd) gl = gl + geom_hline(yintercept = rv$mle_dd[3])
      if(input$CI) gl = gl + geom_errorbar(aes(x=it, y=K, ymin = K-1.96*sdk, ymax = K + 1.96*sdk), colour='darkgreen') 
      gl  + theme_emphasis + geom_hline(yintercept = mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)]), colour="blue")
      
    }
  })
  
  output$fhat <- renderPlot({
    if(length(rv$x)>5){ 
      rv$mcem_it$fhat = rv$fhat
     # if(length(rv$fhat)>9) rv$mcem_it$sdk[10:length(rv$fhat)] = rv$sdk
      gl = ggplot(rv$mcem_it) + geom_line(aes(it,log(fhat))) # + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      if(input$ddd) gl = gl + geom_hline(yintercept = rv$mle_dd[4])
     # if(input$CI) gl = gl + geom_errorbar(aes(x=it, y=K, ymin = K-1.96*sdk, ymax = K + 1.96*sdk), colour='darkgreen') 
      gl  + theme_emphasis + ggtitle(label = "Estimated log likelihood")
      # + geom_hline(yintercept = mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)]), colour="blue")
      
    }
    # if(length(rv$x)>0){
    #   htit <- paste("Estimated loglikelihood, ",length(rv$fhat)," iteratons.")
    #   if(input$log){
    #     plot(input$charts:length(rv$fhat),log(rv$fhat)[input$charts:length(rv$fhat)],col="blue",type="l",main=htit,xlab="EM iteration",ylab="Estimated log-likelihood")
    #     if(input$ddd) points(input$charts:length(rv$ftrue),log(rv$ftrue)[input$charts:length(rv$fhat)])
    #     lines(input$charts:length(rv$fhat),log(rv$fhat+1.96*rv$se)[input$charts:length(rv$fhat)],col="red")
    #     lines(input$charts:length(rv$fhat),log(rv$fhat-1.96*rv$se)[input$charts:length(rv$fhat)],col="red")
    #   }else{
    #     plot(input$charts:length(rv$fhat),rv$fhat[input$charts:length(rv$fhat)],col="blue",type="l",main=htit,xlab="EM iteration",ylab="Estimated likelihood")
    #     if(input$ddd) points(input$charts:length(rv$ftrue),rv$ftrue[input$charts:length(rv$fhat)])
    #     lines(input$charts:length(rv$fhat),(rv$fhat+1.96*rv$se)[input$charts:length(rv$fhat)],col="red")
    #     lines(input$charts:length(rv$fhat),(rv$fhat-1.96*rv$se)[input$charts:length(rv$fhat)],col="red")
    #   }
   # }
    
  })
  

  
  output$diff_fhat_hist <- renderPlot({
    if(length(rv$diff_fhat)>10){ 
      rv$mcem_it$diff_fhat = c(0,rv$diff_fhat)
      gl = ggplot(rv$mcem_it[input$charts:length(rv$diff_fhat),]) + geom_histogram(aes(diff_fhat),binwidth = 0.001)# + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      # if(input$ddd) gl = gl + geom_hline(yintercept = rv$mle_dd[4])
      #gl = ggplot(rv$mcem_it) + geom_histogram(aes(diff_fhat))# + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      
      # if(input$CI) gl = gl + geom_errorbar(aes(x=it, y=K, ymin = K-1.96*sdk, ymax = K + 1.96*sdk), colour='darkgreen') 
      gl  + theme_emphasis + ggtitle(label = "Delta likelihood")
    }
  })
  
  output$rellik_hist <- renderPlot({
    if(length(rv$rellik)>10){ 
    rv$mcem_it$rellik = rv$rellik
    gl = ggplot(rv$mcem_it[input$charts:length(rv$diff_fhat),]) + geom_histogram(aes(rellik),binwidth = 0.0001)# + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
   # if(input$ddd) gl = gl + geom_hline(yintercept = rv$mle_dd[4])
    # if(input$CI) gl = gl + geom_errorbar(aes(x=it, y=K, ymin = K-1.96*sdk, ymax = K + 1.96*sdk), colour='darkgreen') 
    gl  + theme_emphasis + ggtitle(label = "Relative likelihood")
    }
  })
  
  output$hist_w <- renderPlot({
    if(length(rv$x)>0){
      if(input$log_w==TRUE){
        qplot(rv$dim/2,rv$logweights)
      }else{
        qplot(rv$dim/2,rv$weights)
      }
    }
    
  })
}
shinyApp(ui, server)
if (file.exists("first.R")) file.remove("first.R")
data("branching_times_7phylogenies")
n_cores = detectCores()
rv = NULL #Do I need this?
theme_emphasis =  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))  #

MCEM = data.frame(par1=NULL,par2=NULL,par3=NULL,fhat=NULL,fhat.se=NULL,E_time=NULL,M_time=NULL,mc.samplesize=NULL,effective.size=NULL,hessian.inv=NULL)
#pars = c(0,0,0)
#3C0C62  purple
#D8FCFC  lightblue 
##008000 green 
ui <- fluidPage(
  tags$head(tags$style(
    HTML('
         #sidebar {
         background-color:#FDFDFD;
         }
         
         body, label, input, button, select { 
         font-family: "Times New Roman";
         }')
  )),
  sidebarLayout(position = "left",
                
                sidebarPanel(id="sidebar",
                             img(src='emphasis_logo.png', align = "center",height = 170),
                             
                             
                             
                             h3("Data"),  
                             selectInput("brts", "Choose Phylo/Branching times:",
                                         list(
                                           "Dendroica" = "5,4.806886544,4.70731246478,4.50735197578,4.37856240588,4.29594855558,4.19207515688,4.18261061218,4.11238451758,4.09640902445,3.81723693538,3.71143733895,3.48845298905,3.25729503338,3.11613886835,2.64829864145,2.63531839038,2.37990087748,1.82721570435,0.83704715535,0.64242044758,0.56121103655,0.356333544350001,0.346462849050001",
                                           "Anolis" = "103.31057277,97.96837055,94.923127866,94.866216796,90.810203432,90.44410561,90.080885176,86.938065219,83.192481566,79.903508082,78.144291981,75.916079896,75.270062039,74.19732113,72.825735377,72.5234868711,68.360444962,64.1335159681,63.557121926,63.523319671,63.398403586,60.3209181541,59.490443993,57.576137962,56.933279789,56.6964480574,56.361043545,55.506578472,54.495983232,54.199692129,54.051109589,53.692672898,52.897212802,52.23186871,52.164805405,51.94779018,50.553819579,48.373129976,47.4174457904,47.189946167,46.7942740811,44.287638517,43.296282982,43.242616701,42.272773859,42.1266648041,41.905974158,41.329061036,41.1257958974,39.697767108,39.677765636,39.397083778,37.911582502,37.397487349,34.605557178,32.824929892,32.228763421,31.561908554,30.308206955,30.281651159,30.1639183904,29.8042173411,29.773786118,29.6447104204,29.541373926,29.5407793691,28.623740578,28.506256108,27.105138539,26.439039467,26.378922306,26.285718299,26.233564599,24.159222149,22.974342026,21.370573573,21.247374251,20.091077714,19.6384466391,19.483057152,18.985702599,16.272845218,15.9237660074,15.8035460904,15.7840819411,15.2801426021,15.08030346650,13.1506065141,12.572643169,11.2235612480000,10.5640701490000,10.30031894,9.651984572,9.577378633,9.50317724299997,8.43806557499997,6.837926447,5.73566929399999,5.38425275100001,4.4685701345,4.33572126899998,0.828930356499995,0.552543471999996",
                                           "Cetacea" = "35.857845,33.799004,32.390661,31.621529,28.000001,26.063017,26.000001,24.698215,22.044392,19.195664,18.226421,18.023412,17.939427,17.890656,16.066686,15.669702,15.099325,14.540283,14.061555,13.042869,12.847396,11.382959,11.079292,11.028304,10.70277,10.472235,9.438013,8.925942,8.81602,8.803019,8.716039,8.252102,8.20904999999999,8.143058,8.100266,7.677423,7.514394,7.176679,6.975185,6.37563,6.28945,6.04702999999999,5.897506,5.796585,5.616381,5.49323999999999,5.466433,5.27807199999999,5.26532599999999,5.263495,5.096383,4.985876,4.947171,4.927157,4.732447,4.57089299999999,4.45271899999999,4.35571699999999,4.32202299999999,4.170967,4.166225,4.045704,3.791853,3.70627600000000,3.62611499999999,3.44535999999999,3.29116399999999,3.21256099999999,3.07916999999999,3.04867399999999,2.919779,2.83297999999999,2.19441299999999,2.09621900000000,1.93481199999999,1.82123899999999,1.622022,1.570433,1.50673699999999,1.47078099999999,1.36135800000000,1.26811400000000,1.01116300000000,0.924861999999997,0.347030000000004,0.283069999999995",
                                           "Heliconius" = "16.761439,15.160156,14.40605,13.815308,13.486476,13.164424,12.373104,10.840648,10.142988,9.911296,9.273213,9.264573,9.142266,8.536825,8.441098,8.17086,7.92524,7.478269,7.255542,6.851304,5.335066,5.335061,5.152996,4.643518,4.506785,4.446959,3.780976,3.768737,3.488772,3.398945,2.433015,2.048552,1.930075,1.602332,1.302335,1.03376100000000,0.884346",
                                           "Plethodon" = "11.3,9.55365380008198,9.26434040327225,8.83592350352767,8.3434446982257,8.17781491491496,7.95978190384214,6.61207494082374,6.5679856688767,6.21838471981418,5.59809615547134,5.37012669852355,4.7638222125791,4.10749650972075,4.02367324807484,3.65931960175062,3.32916292401100,3.23132222435799,3.18206288699248,2.8572235287017,2.58222342582278,2.43078192215161,1.87377417677032,1.79734091086791,1.77566693721338,1.52675067868777,1.11116172787207,0.800771123741394,0.498973146096477",
                                           "Vangidae" = "9.77,9.652180854,8.612705507,7.491360279,4.94075617,2.5828624",
                                           "Foraminifera" = "64.95,64.9,61.2,44.15,37.2,30.2,23.8,22.5,21,18.8,18.3,17.3,17,14.8,10.8,10.2,10,8.2,8,7.2,5.2,4.5,3.8,3.5,3.4,3.2,2,2,0.8,0.3,0.3",
                                           "Other" = 0
                                           #   "Birds" = paste(as.character(brts_birds), collapse=", "))
                                         )       
                             ),
                             textInput('vec1', 'Or enter a vector (comma delimited) with branching times (Selecting Other)', "4,3.9,3.8,1"),
                             #  actionButton("goTree","Load Tree"),
                             h3("Initial parameters"),
                             numericInput("par1", "Initial lambda_0:", 2),
                             numericInput("par2", "Initial lambda_1:", 0.03),
                             numericInput("par3", "Initial mu_0:", 0.05),
                             
                             h3("Settings"),
                             selectInput("importance_sampler", "Choose Data augmentation sampler:",
                                         list("emphasis" = "emphasis",
                                              "uniform" = "uniform"
                                              )),
                             selectInput("method", "Choose NHPP method:",
                                         list("inverse" = "inverse",
                                              "thinning" = "thinning"
                                         )),
                             selectInput("model", "Choose diversification model:",
                                         list("Diversity dependence" = "dd",
                                              "Exponential diversity dependence" = "edd",
                                              "Phylodiversity dependence" = "pd",
                                              "Exponential Phylodiversity dependence" = "epd",
                                              "GDDX" = "gddx",
                                              "GDPX" = "gpdx"
                                         )),                
                             numericInput("sample_size", "Monte-Carlo sample size:", 10),
                             numericInput("proportion_of_subset", "Proportion of best trees to take:", 1),
                             numericInput("maxspec", "Maximum number of missing species:", 30),
                             h3("Options"),
                        #     conditionalPanel(length(rv$fhat)>10, checkboxInput("CI", "Check CI (after it 10)", FALSE)),
                             #checkboxInput("log", "show estimated lkelihood on log scale", FALSE),
                             checkboxInput("log_w", "show weights on log scale", FALSE),
                             numericInput("charts", "See charts from iteration:", 1),
                       
                        #actionButton("resetbutt","Reset"),
                             checkboxInput("parallel", "Parallel", TRUE),
                             numericInput("cores",paste("Your computer holds",n_cores,"cores, how many of them you want to use?"),2),
                       h3("Controls"),
                       actionButton("gogobutt","Go"),
                       actionButton("stopbutt","Stop"),
                       h3("Save"),
                       checkboxInput("save", "Save Current MCEM state", FALSE),
                       textInput("file", "Directory where you want to save:", "~/Google Drive/scripts for jobs and data/Experiments/MCEM/MCEM.RData")
                       
                ),
                mainPanel(
                  tabsetPanel(type = "tabs",
                              
                              tabPanel("Help",
                                       h1("Welcome to Emphasis!"),
                                       
                                     #  "Emphasis mework for a diversity dependence model.\n",
                                       "This UI version of the emphasis framework allows interactive and user friendly usage for Emphasis analysis.\n",
                                       "If you are familiar with the emphasis approach you just need a short description of the functionalities of this app and we are good to go.\n",
                                     #  "As a simple example we suggest to type the tree 4,3.9,3.8,1,0.5.",
                                       
                                       h3("Data"),
                                       "Load your phylogenetic tree. At the moment (for current models) data consists on a vector with branching times starting from stem age in decreasing order.",
                                       "It is also possible to use one of the 7 data examples",
                                       "Below we vizualize the tree to analyze",
                                       plotOutput("phylogenetic_tree"),
                                       
                                      
                                      h3("Diversification model"),
                                      "Select the model for diversification rates.",
                                      
                                      # textOutput("diversification_model"),
                                    
                                       includeHTML("dd.html"),
                                     
                                     
                                       h3("Initial parameters"),
                                       "Initial parameters for diversity dependence model.",
                                       h3("Settings"),
                                       h4("Data augmentation sampler"),
                                       "Importance sampler to be used on the data augmentation scheme.",
                                       h4("Diversification model"),
                                       "Diversification rates model. Currently constan-rates and diversity-dependence are available.",
                                       h4("Monte Carlo Sampling size"),
                                       "This is the Monte Carlo sampling size on the E step. You can check if it is large enough looking at the variation on the plots as well as the weights plot. It is possible to increase it at any moment of the routine, the change will apply after the current iteration next iteration.",
                                       h4("Number of best trees to take"),
                                       "Numbers of trees to take on the M-step. Check with the weights plot and the percentage of likelihood information if the current is appropriate or not.",
                                       h4("Maximum number of species"),
                                       "The current data augmentation algorithm requires a limit on the number of missing species, you can check if that limit makes sense by checking the weights plot. If the weights are spread all over the plot you probably need to increase the maximum number of species allowed to sample in the whole space. If the last meaningful weight is far to the left to your limit then you probably want to decrease this quantity for the sake of efficiency.",
                                       h3("Options"),
                                       "Several options are available for analysis",
                                       #  h4("Compare with DDD"),
                                       #   "For the given examples (except birds\' tree) a solution is given by DDD package, useful for comparison.",
                                       #  h4("log of estimated likelihood"),
                                       #"To observe the estimated loglikelihood",
                                       h4("log of weights"),
                                       "To observe the weights on logarithm scale",
                                       #h4("Iteration to show on plot"),
                                      # "To observe the plots on more detail from last iterations",
                                     h3("Controls"),
                                     "To initialize and stop the MCEM routine.",
                                       h4("Number of processors"),
                                       "Number of cores to be used on parallel. We suggest 2 for small samples and increase it as the sample size increases. Ready to go to the analysis tab?"
                              ),
                              
                              tabPanel("Analysis",
                                       
                                       fluidRow(
                                         #dataTableOutput("iterations_data"),
                                     #    paste(c("current parameters: ",pars)),
                                        
                                         column(5,
                                    #    
                                                h3("Parameters"),
                                                plotOutput("lambda"),
                                                plotOutput("mu"),
                                                plotOutput("K")),
                                         
                                  #       column(5,
                                   #             h3("Diagnostics"),
                                   #             plotOutput("fhat"),
                                   #             plotOutput("diff_fhat_hist")#,
                                   #             # plotOutput("rellik_hist")
                                   #      ),
                                         column(5,
                                                h3("Weights"),
                                                plotOutput("weights_by_dimension"),  #weight_vs_dimension
                                                plotOutput("fvsg")
                                         ),
                                         column(5,
                                                h3("Information"),
                                                plotOutput("timeconsumption"),
                                              #  textOutput("txtOutput2"),
                                                textOutput("txtOutput3")
                                         )
                                       )),
                              
                              tabPanel("About",
                                       h2(""),
                                       #img(src='logo.png', align = "left",height = 180),
                                       helpText("Emphasis is a joint work by Francisco Richter, Rampal Etienne and Ernst Wit. For questions or comments please write to f.richter@rug.nl.\n", align = "left",height = 70),
                                       img(src='license.png', align = "right",height = 60)
                              )
                  ) #close tab panel
                ) #close mainpanel
                
  ) # close sidebarLayout
) #close UI






server <- shinyServer(function(input,output,session) {
  
  # output$open_end <- renderUI({
  #   if (!(input$brts[1] == 0)){
  #     return(NULL)
  #   }else {
  #     textInput('vec1', 'Or enter a vector (comma delimited) with branching times (Selecting Other)', "4,3.9,3.8,1")
  #   }
  #  })
  
  
  
  rv <- reactiveValues(MCEM=MCEM,
                       run=F,
                       sdl=NULL,
                       sdm=NULL,
                       sdk=NULL,
                       dim=NULL,
                       dims=NULL,
                       weights=NULL,
                       weights_by_dims=NULL,
                       logweights=NULL,
                       em.iteration=0,
                       logf=NULL,
                       logg=NULL)
  input_values <- reactiveValues(brts=NULL,init_pars=NULL,brts_d=NULL)
  autoInvalidate <- reactiveTimer(intervalMs=1000,session)
  
  observe({
    # load branching times
    if(input$brts[1]==0){
      brts = brts_d <- as.numeric(unlist(strsplit(input$vec1,",")))
    }else{ 
      brts = brts_d <- as.numeric(unlist(strsplit(input$brts,",")))
    }
    # load reactive input values
    input_values$brts = brts
    input_values$init_pars = c(input$par1,input$par2,input$par3)
    
  })
  
  observe({
    pars = input_values$init_pars
    autoInvalidate()  # to make it run withouth changing input 
    isolate({ if (rv$run) { 
      
      if(file.exists("first.R")) load("first.R")
      rv$em.iteration = rv$em.iteration + 1
      
      withProgress(message = paste("Performing MCEM: iteration ", rv$em.iteration), value = 0, {
        
        time = proc.time()
        setProgress(value=0,detail = "Performing E step")
        # st = mc_augmentation_thinning(brts=input_values$brts,pars = pars,model = input$model,importance_sampler = input$method,sample_size = input$sample_size,parallel = FALSE,no_cores = input$cores)
        st = mc_sample_independent_trees(brts = input_values$brts,pars = pars,nsim = input$sample_size,model = input$model, importance_sampler = input$importance_sampler,no_cores = input$cores,maxnumspec = input$maxspec,method=input$method)
        E_time = get.time(time)
        
        time = proc.time()
        setProgress(value=0.5,detail = "Performing M step")
        M = M_step(st = st,init_par = pars,model = input$model,proportion_of_subset=input$proportion_of_subset)
        M_time = get.time(time)
        # if(!is.na(M$po$value)){
          pars = M$po$par
        # }
        # hessian_inverse = try(diag(solve(M$po$hessian)))
        fhat = st$fhat
        # se = st$fhat.se
        # if(!is.numeric(h1)) h1 = c(NULL,NULL,NULL)
        
        mcem = list(pars=pars,fhat=fhat,st=st,effective_sample_size=M$effective_sample_size,E_time=E_time,M_time=M_time)
        
      })
      
  
      MCEM_last = data.frame(par1=pars[1],
                             par2=pars[2],
                             par3=pars[3],
                             fhat=mcem$fhat,
                           #  fhat.se=mcem$se,
                             E_time=mcem$E_time,
                             M_time=mcem$M_time,
                             mc.samplesize=input$sample_size,
                             cores = input$cores,
                             effective.size=mcem$effective_sample_size,
                           #  hessian.inv1 = mcem$hessian_inverse[1],
                          #   hessian.inv2 = mcem$hessian_inverse[2],
                          #   hessian.inv3 = mcem$hessian_inverse[3],
                             sdl = NaN,
                             sdm = NaN,
                             sdk = NaN,
                             em.iteration = rv$em.iteration,
                             diff_fhat = NaN
      )
      rv$MCEM = rbind(rv$MCEM,MCEM_last)
      
      rv$weights = mcem$st$weights
      rv$logweights = mcem$st$logweights
      rv$dim = sapply(mcem$st$trees,FUN = function(list) nrow(list))
      
      ###
      ta = table(rv$dim)
      dims = as.numeric(names(ta))
      lenghts = as.numeric(ta)
      weights = NULL
      for(i in 1:length(dims)){
        weights[i] = sum(mcem$st$weights[rv$dim == dims[i]])
      }
      rv$weights_by_dims=weights/input$sample_size
      rv$dims=dims
      ###
      rv$logf = mcem$st$logf
      rv$logg = mcem$st$logg
      previous_pars=pars
      pars = mcem$pars
      
      save(pars,previous_pars,file="first.R")
      MCEM_temp = rv$MCEM
      
      if(input$save){
        withProgress(message = 'Saving current state', value = 0, {
          save(MCEM_temp,file=input$file)
          incProgress(3/10, detail = paste("Saving Data", 2))
        })
      }  
      
    } })#end isolate
  }) #end observe
  
  observeEvent(input$gogobutt, { isolate({ rv$run=T      }) })
  observeEvent(input$stopbutt, { isolate({ rv$run=F      }) })
  
  
  output$diversification_model = renderText({
    if(input$model=="dd"){
      includeMarkdown(rmarkdown::render("dd.Rmd"))
    }
    if(input$model=="edd"){ 
      includeHTML("edd.html")
    }
    #includeMarkdown("dd.Rmd")
    includeHTML("edd.html")
  })
  
  output$iterations_data <- DT::renderDataTable(DT::datatable({
    rv$MCEM
  }))
  

  
  output$txtOutput3 = renderText({
    paste0("Proportion of likelihood: ", rv$MCEM$effective.size[nrow(rv$MCEM)] )
  })
  
  

  output$lambda <- renderPlot({
    
    if(nrow(rv$MCEM)>2){ 
      lambda.est = rv$MCEM$par1[nrow(rv$MCEM)]
      MCEM = rv$MCEM[input$charts:nrow(rv$MCEM),]
      plot.lambda = ggplot(MCEM) + geom_line(aes(em.iteration,par1)) + ggtitle(label=paste("Last estimation:  ", lambda.est)) + ylab("parameter 1") + xlab("EM iteration")
    #  if(input$CI & nrow(rv$MCEM)>12) plot.lambda = plot.lambda + geom_errorbar(aes(x=em.iteration, y=par1, ymin = par1-1.96*sdl, ymax = par1 + 1.96*sdl), colour='darkgreen') 
      change.ss = which(diff(rv$MCEM$mc.samplesize)>0)
      if(sum(change.ss)>0){
        plot.lambda = plot.lambda + geom_vline(xintercept = change.ss, colour="darkgreen",linetype="dotted")
      }
      plot.lambda + theme_emphasis
    }
  })
  
  output$mu <- renderPlot({
    if(nrow(rv$MCEM)>2){ 
      mu.est = rv$MCEM$par2[length(rv$MCEM$par2)]
      MCEM = rv$MCEM[input$charts:length(rv$MCEM$par3),]
      plot.mu = ggplot(MCEM) + geom_line(aes(em.iteration,par2)) + ggtitle(label=paste("Last estimation:  ",mu.est)) + ylab("Parameter 2") + xlab("EM iteration")
    #  if(input$CI & nrow(rv$MCEM)>12) plot.mu = plot.mu + geom_errorbar(aes(x=em.iteration, y=par2, ymin = par2-1.96*sdm, ymax = par2 + 1.96*sdm), colour='darkgreen') 
      change.ss = which(diff(rv$MCEM$mc.samplesize)>0)
      if(sum(change.ss)>0){
        plot.mu = plot.mu + geom_vline(xintercept = change.ss, colour="darkgreen",linetype="dotted")
      }
      plot.mu + theme_emphasis
    }
  })
  
  output$K <- renderPlot({
    if(nrow(rv$MCEM)>2){ 
      K.est = mean(rv$MCEM$par3[input$charts:length(rv$MCEM$par3)])
      MCEM = rv$MCEM[input$charts:length(rv$MCEM$par3),]
      K.plot = ggplot(MCEM) + geom_line(aes(em.iteration,par3)) + ggtitle(label=paste("Last estimation:  ",K.est)) + ylab("Parameter 3") + xlab("EM iteration")
      change.ss = which(diff(rv$MCEM$mc.samplesize)>0)
      if(sum(change.ss)>0){
        K.plot = K.plot + geom_vline(xintercept = change.ss, colour="darkgreen",linetype="dotted")
      }
      K.plot  + theme_emphasis 
    }
  })
  
  output$fhat <- renderPlot({
    if(nrow(rv$MCEM)>5){ 
      MCEM = rv$MCEM[input$charts:length(rv$MCEM$par3),]
      fhat.plot = ggplot(MCEM) + geom_line(aes(em.iteration,log(fhat))) # + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      #  if(input$ddd) fhat.plot = fhat.plot + geom_point(aes(em.iteration,ftrue))
      if(input$CI) fhat.plot = fhat.plot + geom_line(aes(em.iteration,log(rv$MCEM$fhat+1.96*rv$MCEM$fhat.se)),col="red") + geom_line(aes(em.iteration,log(rv$MCEM$fhat-1.96*rv$MCEM$fhat.se)),col="red")#+ geom_errorbar(aes(x=it, y=K, ymin = K-1.96*sdk, ymax = K + 1.96*sdk), colour='darkgreen')
      change.ss = which(diff(rv$MCEM$mc.samplesize)>0)
      if(sum(change.ss)>0){
        fhat.plot = fhat.plot + geom_vline(xintercept = change.ss, colour="darkgreen",linetype="dotted")
      }
      fhat.plot  + theme_emphasis + ggtitle(label = "Estimated log likelihood")
    }
  })
  
  output$diff_fhat_hist <- renderPlot({
    if(nrow(rv$MCEM)>12){ 
      breaks <- pretty(range(rv$MCEM$diff_fhat), n = nclass.FD(rv$MCEM$diff_fhat), min.n = 1)
      bwidth <- breaks[2]-breaks[1]
      rv$MCEM$diff_fhat = c(0,diff(rv$MCEM$fhat))
      gl = ggplot(rv$MCEM[input$charts:nrow(rv$MCEM),]) + geom_histogram(aes(diff_fhat),binwidth = bwidth)#,binwidth = (max(diff_fhat)-max(diff_fhat))/50)# + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      gl  + theme_emphasis + ggtitle(label = "Delta likelihood")
      #stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1))
    }
  })
  
  output$timeconsumption <- renderPlot({
    if(nrow(rv$MCEM)>2){ 
      df = data.frame(it=rep(rv$MCEM$em.iteration,2),time=c(rv$MCEM$E_time,rv$MCEM$M_time),Process=c(rep("E_step",nrow(rv$MCEM)),rep("M_step",nrow(rv$MCEM))))
      time_consumption_plot = ggplot(df, aes(x=it, y=time, fill=Process)) + geom_area()
      change_up = which(diff(rv$MCEM$cores)>0)
      change_down = which(diff(rv$MCEM$cores)<0)
      if(sum(change_up)>0){
        time_consumption_plot = time_consumption_plot + geom_vline(xintercept = change_up, colour="darkgreen",linetype="dotted")
      }
      if(sum(change_down)>0){
        time_consumption_plot = time_consumption_plot + geom_vline(xintercept = change_down, colour="green",linetype="dotted")
      }
      time_consumption_plot + ggtitle("Time consumption per MCEM iteration") + theme_emphasis 
    }
  })
  
  output$weights_by_dimension <- renderPlot({
    if(nrow(rv$MCEM)>2){
      DF = data.frame(dims=rv$dims,weights=rv$weights_by_dims)
      if(input$log_w){
        weights_by_dimension = ggplot(DF) + geom_line(aes(x=dims,y=log(weights)))
      }else{
        weights_by_dimension = ggplot(DF) + geom_line(aes(x=dims,y=weights))
      }
      weights_by_dimension + theme_emphasis 
    }
  })
  
  output$phylogenetic_tree <- renderPlot({
    tree = list(wt = -diff(c(input_values$brts,0)),to=rep(1,length(input_values$brts)-1))
    tree = vectors2phylo(tree)
    plot(tree,show.tip.label=F,root.edge=TRUE)
    axisPhylo()
  })
  
  output$fvsg <- renderPlot({
    if(nrow(rv$MCEM)>2){
      datita = data.frame(logf=rv$logf,logg=rv$logg,dim=rv$dim)
      ggplot(datita) + geom_point(aes(x=logf,y=logg,colour=dim)) + theme_emphasis
    }
  })
}
)

shinyApp(ui, server)
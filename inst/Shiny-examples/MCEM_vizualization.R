library(shiny)
library(ggplot2)
rm(list = ls())


# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Emphasis Analytics"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      fileInput("file1", "Choose RData File",
                multiple = TRUE,
                accept = ".RData"
      ),
      
      # Input: Selector for choosing dataset ----
      selectInput(inputId = "par",
                  label = "Choose plot parameter:",
                  choices = c("par1", "par2", "par3","par1VSpar3","eff.sizeVSpar1","oneline","lambdaHist")),
      
      uiOutput('columns1'),
      uiOutput('columns2'),
      
      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "obs",
                   label = "Number of observations to view:",
                   value = 10),
      
      numericInput(inputId = "burning",
                   label = "Number of observations to burn:",
                   value = 10),
      
      numericInput(inputId = "filter",
                   label = "Filter by size:",
                   value = 0),
      
      selectInput(inputId = "rep",
                  label = "Choose plot replicant:",
                  choices = c("all",1:2))
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      plotOutput("parameter_estimation_general"),
      
      # Output: Verbatim text for data summary ----
      verbatimTextOutput("summary")#,
      #verbatimTextOutput("summary2"),
     # verbatimTextOutput("In")
      
      # Output: HTML table with requested number of observations ----
      #tableOutput("view")
      #plotOutput("sample_size"),
      #plotOutput("lambda2"),
      #plotOutput("logfhat")
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  # Return the requested dataset ----
  # By declaring datasetInput as a reactive expression we ensure
  # that:
  #
  # 1. It is only called when the inputs it depends on changes
  # 2. The computation and result are shared by all the callers,
  #    i.e. it only executes a single time
 
  
  DF <- reactive({
    inFile <- input$file1
    DF_temp = NULL
    for(i in 1:nrow(inFile)){
      load(inFile[[i, 'datapath']])
      if(exists("MCEM_temp")){
        DF = data.frame(iteration=MCEM_temp$em.iteration,par1=MCEM_temp$par1,par2=MCEM_temp$par2,par3=MCEM_temp$par3,E_time=MCEM_temp$E_time,Mtime=MCEM_temp$M_time,efective_sample_size=MCEM_temp$effective.size,h1=1/MCEM_temp$hessian.inv1,h2=1/MCEM_temp$hessian.inv2,h3=1/MCEM_temp$hessian.inv3)
        rm(MCEM_temp)
      }
      if("pars3" %in% names(DF)){
        names(DF)[names(DF)=="pars3"] <- "par3"
      }
      DF$rep = i
      DF_temp = rbind(DF_temp,DF)
    }
    DF_temp$rep = as.character(DF_temp$rep)
    return(DF_temp)
  })
  
  In <- reactive({
    inFile <- input$file1
    load(inFile$datapath)
    return(input)
  })
  
  output$columns1 = renderUI({
    selectInput(inputId = "par1",
                label = "Choose x plot parameter:",
                choices=names(DF()))
  })
  
  output$columns2 = renderUI({
    selectInput(inputId = "par2",
                label = "Choose y plot parameter:",
                choices=names(DF()))
  })
  
  
  
  # Create caption ----
  # The output$caption is computed based on a reactive expression
  # that returns input$caption. When the user changes the
  # "caption" field:
  #
  # 1. This function is automatically called to recompute the output
  # 2. New caption is pushed back to the browser for re-display
  #
  # Note that because the data-oriented reactive expressions
  # below don't depend on input$caption, those expressions are
  # NOT called when input$caption changes
  output$caption <- renderText({
    input$caption
  })
  
  # Generate a summary of the dataset ----
  # The output$summary depends on the datasetInput reactive
  # expression, so will be re-executed whenever datasetInput is
  # invalidated, i.e. whenever the input$dataset changes
  output$summary <- renderPrint({
    df = DF()
    df = df[input$burning:nrow(df),]
    summary(df)
  })
  
  output$summary2 <- renderPrint({
    df = DF()
    df = df[input$burning:nrow(df),]
    df = df[df$efective_sample_size > input$filter,]
    summary(log(df$fhat))
  })
  
  output$In <- renderPrint({
    In()
  })
  
  # Show the first "n" observations ----
  # The output$view depends on both the databaseInput reactive
  # expression and input$obs, so it will be re-executed whenever
  # input$dataset or input$obs is changed
  output$view <- renderTable({
    head(datasetInput(), n = input$obs)
  })
  
  output$parameter_estimation <- renderPlot({
    df = DF()
    df = df[input$obs:nrow(df),]
    df = df[df$efective_sample_size >= input$filter,]
    if(input$par == "par1"){
      par_est = mean(df$par1[input$burning:nrow(df)])
      plot_par_est = ggplot(df)+geom_line(aes(x=iteration,y=par1, colour=rep))+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$par == "par2"){
      par_est = mean(df$par2[input$burning:nrow(df)])
      plot_par_est = ggplot(df)+geom_line(aes(x=iteration,y=par2))+ggtitle("beta Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$par == "par3"){
      par_est = mean(df$par3[input$burning:nrow(df)])
      plot_par_est = ggplot(df)+geom_line(aes(x=iteration,y=par3))+ggtitle("mu Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$par == "par1VSpar3"){
      plot_par_est = ggplot(df)+geom_path(aes(x=par1,y=par3,colour=rep))+ggtitle("lambda vs mu Estimation")#+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$par == "eff.sizeVSpar1"){
      plot_par_est = ggplot(df)+geom_path(aes(x=efective_sample_size,y=par1,colour=rep))+ggtitle("lambda vs mu Estimation")#+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$par == "oneline"){
      df2 = df[order(df$efective_sample_size),]
      df2$iteration = 1:nrow(df2)
      plot_par_est = ggplot(df2) + geom_line(aes(x=iteration,y=par1))
    }
    if(input$par == "lambdaHist"){
      plot_par_est = ggplot(df)+geom_histogram(aes(x=par1))
    }
    if(input$par == "samplesize"){
      plot_par_est = ggplot(df)+geom_line(aes(x=iteration,y=efective_sample_size))+ggtitle("Effective sample size")
    }
    
    
    
    if("mc.samplesize" %in% names(df)){
      change.ss = which(diff(df$mc.samplesize)>0)
    }else{
      change.ss = NULL
    }
    if(sum(change.ss)>0){
      plot_par_est = plot_par_est + geom_vline(xintercept = change.ss, colour="darkgreen",linetype="dotted")
    }
    plot_par_est
  })
  
  output$parameter_estimation_general <- renderPlot({
    df = DF()
    df = df[input$obs:nrow(df),]
    df = df[df$efective_sample_size >= input$filter,]
   # df2 = data.frame()
    plot_par_est = ggplot(df,aes(colour=rep))+geom_point(aes_string(x=input$par1,y=input$par2))#+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    plot_par_est
  })
  
  
  output$sample_size <-renderPlot({
    df = DF()
    df = df[input$obs:nrow(df),]
    change.ss = which(diff(df$mc.samplesize)>0)
    plot.lambda = ggplot(df)+geom_line(aes(x=iteration,y=efective_sample_size))+ggtitle("Effective sample size")#+geom_hline(yintercept = lambda.est,colour="red")
    if(sum(change.ss)>0){
      plot.lambda = plot.lambda + geom_vline(xintercept = change.ss, colour="darkgreen",linetype="dotted")
    }
    plot.lambda
  })
  
  output$lambda2 <- renderPlot({
    df = DF()
    lambda.est = mean(df$par1[input$burning:nrow(df)])
    df = df[input$obs:nrow(df),]
    ggplot(df)+geom_histogram(aes(x=par1))

  })
  
  output$logfhat <- renderPlot({
    df = DF()
    #lambda.est = mean(df$par1[input$burning:nrow(df)])
    df = df[input$obs:nrow(df),]
    ggplot(df)+geom_histogram(aes(x=log(fhat)))
    
  })
  
  output$mu <- renderPlot({
    if(nrow(rv$MCEM)>5){ 
      mu.est = rv$MCEM$par2[length(rv$MCEM$par2)]
      MCEM = rv$MCEM[input$charts:length(rv$MCEM$par3),]
      plot.mu = ggplot(MCEM) + geom_line(aes(iteration,par2)) + ggtitle(label=paste("Last estimation:  ",mu.est)) + ylab(expression(lambda1)) + xlab("EM iteration")
      if(input$CI & nrow(rv$MCEM)>12) plot.mu = plot.mu + geom_errorbar(aes(x=iteration, y=par2, ymin = par2-1.96*sdm, ymax = par2 + 1.96*sdm), colour='darkgreen') 
      change.ss = which(diff(rv$MCEM$mc.samplesize)>0)
      if(sum(change.ss)>0){
        plot.mu = plot.mu + geom_vline(xintercept = change.ss, colour="darkgreen",linetype="dotted")
      }
      plot.mu + theme_emphasis
    }
  })
  
  output$K <- renderPlot({
    if(nrow(rv$MCEM)>5){ 
      K.est = mean(rv$MCEM$par3[input$charts:length(rv$MCEM$par3)])
      MCEM = rv$MCEM[input$charts:length(rv$MCEM$par3),]
      K.plot = ggplot(MCEM) + geom_line(aes(iteration,par3)) + ggtitle(label=paste("Last estimation:  ",K.est)) + ylab(expression(mu)) + xlab("EM iteration")
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
      fhat.plot = ggplot(MCEM) + geom_line(aes(iteration,log(fhat))) # + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      #  if(input$ddd) fhat.plot = fhat.plot + geom_point(aes(em.iteration,ftrue))
      if(input$CI) fhat.plot = fhat.plot + geom_line(aes(iteration,log(rv$MCEM$fhat+1.96*rv$MCEM$fhat.se)),col="red") + geom_line(aes(iteration,log(rv$MCEM$fhat-1.96*rv$MCEM$fhat.se)),col="red")#+ geom_errorbar(aes(x=it, y=K, ymin = K-1.96*sdk, ymax = K + 1.96*sdk), colour='darkgreen')
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
  
}

# Create Shiny app ----
shinyApp(ui, server)
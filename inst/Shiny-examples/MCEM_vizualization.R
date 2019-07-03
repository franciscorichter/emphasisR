library(shiny)

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Emphasis Analytics"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      fileInput("file1", "Choose RData File",
                accept = ".RData"
      ),
      
      # Input: Selector for choosing dataset ----
      selectInput(inputId = "dataset",
                  label = "Choose parameter to see:",
                  choices = c("par1", "par2", "par3")),
      
      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "obs",
                   label = "Number of observations to view:",
                   value = 10),
      
      numericInput(inputId = "burning",
                   label = "Number of observations to burn:",
                   value = 10)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      plotOutput("lambda"),
      
      # Output: Verbatim text for data summary ----
      verbatimTextOutput("summary"),
      verbatimTextOutput("summary2"),
      
      # Output: HTML table with requested number of observations ----
      #tableOutput("view")
      plotOutput("lambda2"),
      plotOutput("logfhat")
      
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
    load(inFile$datapath)
    if(exists("MCEM_temp")){
      DF = MCEM_temp
    }else{
      DF$em.iteration = DF$iteration
    }
    return(DF)
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
    summary(log(df$fhat))
  })
  
  # Show the first "n" observations ----
  # The output$view depends on both the databaseInput reactive
  # expression and input$obs, so it will be re-executed whenever
  # input$dataset or input$obs is changed
  output$view <- renderTable({
    head(datasetInput(), n = input$obs)
  })
  
  
  output$lambda <-renderPlot({
    df = DF()
    lambda.est = mean(df$par1[input$burning:nrow(df)])
    df = df[input$obs:nrow(df),]
    change.ss = which(diff(df$mc.samplesize)>0)
    plot.lambda = ggplot(df)+geom_line(aes(x=em.iteration,y=par1))+ggtitle("Lambda Estimation",subtitle = paste("Current estimation:",lambda.est))+geom_hline(yintercept = lambda.est,colour="red")
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
      plot.mu = ggplot(MCEM) + geom_line(aes(em.iteration,par2)) + ggtitle(label=paste("Last estimation:  ",mu.est)) + ylab(expression(lambda1)) + xlab("EM iteration")
      if(input$CI & nrow(rv$MCEM)>12) plot.mu = plot.mu + geom_errorbar(aes(x=em.iteration, y=par2, ymin = par2-1.96*sdm, ymax = par2 + 1.96*sdm), colour='darkgreen') 
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
      K.plot = ggplot(MCEM) + geom_line(aes(em.iteration,par3)) + ggtitle(label=paste("Last estimation:  ",K.est)) + ylab(expression(mu)) + xlab("EM iteration")
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
  
}

# Create Shiny app ----
shinyApp(ui, server)
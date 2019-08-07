library(emphasis)
library(xtable)
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
      uiOutput('columns1'),
      uiOutput('columns2'),
      
      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "obs",
                   label = "Number of plot observations:",
                   value = 1),
      
      selectInput(inputId = "typePlot",
                  label = "Type of plot:",
                  choices = c("Path","Points")),
      
      numericInput(inputId = "burning",
                   label = "Number of observations to burn:",
                   value = 1),
      
      numericInput(inputId = "filter",
                   label = "Set maximum number of iterations to consider for all data sets:",
                   value = 2000),
      
      checkboxInput("logY",
                    label="Log scale on Y axis"),
      
      selectInput(inputId = "rep",
                  label = "Choose plot replicant:",
                  choices = c("all",1:2))
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      plotOutput("parameter_estimation_general"),
      
      # Output: Verbatim text for data summary ----
      
      tableOutput("view"),
      verbatimTextOutput("summary"),
      verbatimTextOutput("code")
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
 
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
    save(DF_temp,file="palPaper.RData")
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
  
  output$caption <- renderText({
    input$caption
  })

  output$summary <- renderPrint({
    df = DF()
    df = df[input$burning:input$filter,]
    summary(df)
  })
  
  data_to_table <- function(df,replicant,init_plot=input$burning,end_point=input$filter){
    df = df[df$rep==replicant,]
    df = df[init_plot:end_point,]
    gampar1 = gam(par1 ~ s(iteration), data=df)
    standard_error_par1 = sqrt(1/df$h1+gampar1$sig2/nrow(df))
    gampar2 = gam(par2 ~ s(iteration), data=df)
    standard_error_par2 = sqrt(1/df$h2+gampar2$sig2/nrow(df))
    gampar3 = gam(par3 ~ s(iteration), data=df)
    standard_error_par3 = sqrt(1/df$h3+gampar3$sig2/nrow(df))
    summ = data.frame(MSS = mean(df$efective_sample_size),replicant=replicant,par1=mean(df$par1),SD1=mean(standard_error_par1),par2=mean(df$par2),SD2=mean(standard_error_par2),par3=mean(df$par3),SD3=mean(standard_error_par3))
    return(summ)
  }
  
  output$code <- renderPrint({
    ta=NULL
    df = DF()
    for(i in unique(df$rep)){
      ta = rbind(ta,data_to_table(df=df,replicant=i,init_plot=input$burning)) 
    }
    xtable(ta,digits=3)
  })
  
  output$In <- renderPrint({
    In()
  })
  
  output$view <- renderTable({
    ta=NULL
    df = DF()
    for(i in unique(df$rep)){
      ta = rbind(ta,data_to_table(df=df,replicant=i,init_plot=input$burning)) 
    }
    ta
  },
   digits = 6)
  
 
  
  output$parameter_estimation_general <- renderPlot({
    df = DF()
    df = df[input$obs:nrow(df),]
    if(input$typePlot=="Path"){
      plot_par_est = ggplot(df,aes(colour=rep))+geom_path(aes_string(x=input$par1,y=input$par2))#+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$typePlot=="Points"){
      plot_par_est = ggplot(df,aes(colour=rep))+geom_point(aes_string(x=input$par1,y=input$par2,alpha=(1:nrow(df))/(2*nrow(df))))#+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$logY){
      plot_par_est = plot_par_est + scale_y_log10() 
    }
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
  
  
  output$logfhat <- renderPlot({
    df = DF()
    #lambda.est = mean(df$par1[input$burning:nrow(df)])
    df = df[input$obs:nrow(df),]
    ggplot(df)+geom_histogram(aes(x=log(fhat)))
    
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
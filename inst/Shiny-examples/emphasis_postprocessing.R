library(emphasis)
library(GGally)
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
      
      
      
      radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                   inline = TRUE),
      downloadButton('downloadReport'),
      # downloadButton("report", "Generate report"),
      
      selectInput(inputId = "typePlot",
                  label = "Type of plot:",
                  choices = c("Path","Points")),
      
      numericInput(inputId = "burn1",
                   label = "Number of observations to burn:",
                   value = 1),
      
      numericInput(inputId = "burn2b",
                   label = "Set maximum number of iterations to consider for all data sets:",
                   value = 200),
      
      checkboxInput("logY",
                    label="Log scale on Y axis"),
      
      selectInput(inputId = "rep",
                  label = "Choose plot replicant:",
                  choices = c("all",1:2))
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      numericInput(inputId = "obs1",
                   label = "left limit for plot:",
                   value = 1),
      numericInput(inputId = "obs2",
                   label = "right limit for plot:",
                   value = 100),
      plotOutput("parameter_estimation_general", click = "plot_click"),
      uiOutput("correlations_tab"),
      plotOutput("correlations"),
      # Output: Verbatim text for data summary ----
      textOutput("caption"),
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
      
      file = inFile[[i, 'datapath']]
      load(file)
      
      
      DF1 = data.frame(iteration=1:nrow(mcem),mu0=mcem$par1,lambda0=mcem$par2,N=mcem$par3,PD=mcem$par4,E_time=mcem$E_time,M_time=mcem$M_time,sample_size=mcem$sample_size,fhat=mcem$fhat,rep=as.character(i))
      
      DF_temp = rbind(DF_temp,DF1)
      
    }
    #save(DF_temp,file="palPaper.RData")
    return(DF_temp)
  })
  
  In <- reactive({
    inFile <- input$file1
    load(inFile$datapath)
    return(input)
  })
  
  right = reactive({
    if(is.null(input$plot_click$x)){
      val = 100
    }else{
      val = input$plot_click$x
    }
    val
  })
  
  output$columns1 = renderUI({
    selectInput(inputId = "par1",
                label = "Choose x plot parameter:",
                choices=names(DF()))
  })
  
  output$columns2 = renderUI({
    selectInput(inputId = "par2",
                label = "Choose y plot parameter:",
                choices=rev(names(DF())))
  })
  
  output$correlations_tab = renderUI({
    selectInput(inputId = "replic",
                label = "Choose replicant:",
                choices = unique(DF()$rep))
  })
  
  output$caption <- renderPrint({
    print(In())
  })
  
  output$summary <- renderPrint({
    df = DF()
    df = df[df$iteration %in% input$obs1:input$obs2,]
    summary(df)
  })
  
  data_to_table <- function(df,replicant,init_plot=input$burn1,end_point=right()){
    df = df[df$rep==replicant,]
    df = df[df$iteration %in% input$obs1:input$obs2,]
    summ = data.frame(lfhat = mean(df$fhat),replicant=replicant,par1=median(df$par1),par2=median(df$par2),par3=median(df$par3),par4=median(df$par4),E_time = median(df$E_time), M_time = median(df$M_time), sample_size=mean(df$sample_size),cores=mean(df$cores))
    return(summ)
  }
  
  output$code <- renderPrint({
    ta=NULL
    df = DF()
    for(i in unique(df$rep)){
      ta = rbind(ta,data_to_table(df=df,replicant=i,init_plot=input$burn1)) 
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
      ta = rbind(ta,data_to_table(df=df,replicant=i,init_plot=input$burn1)) 
    }
    ta
  },
  digits = 6)
  
  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('my-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    },
    
    content = function(file) {
      src <- normalizePath('report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      
      library(rmarkdown)
      out <- render('report.Rmd', switch(
        input$format,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
  )
  
  output$parameter_estimation_general <- renderPlot({
    df = DF()
    if(input$obs2>nrow(df)){
      lim2 = nrow(df)
    }else{
      lim2 = input$obs2
    }
    #df = df[df$iteration %in% input$obs1:lim2 & !is.na(df$mu0),]
    if(input$typePlot=="Path"){
      plot_par_est = ggplot(df,aes(colour=rep))+geom_path(aes_string(x=input$par1,y=input$par2))#+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$typePlot=="Points"){
      plot_par_est = ggplot(df,aes(colour=rep))+geom_point(aes_string(x=input$par1,y=input$par2,alpha=(1:nrow(df))/(2*nrow(df))))#+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
 #   if(input$logY){
#      plot_par_est = plot_par_est + scale_y_log10() 
#    }
    if(input$par1=="iteration"){
      plot_par_est = plot_par_est +geom_vline(xintercept = input$obs1)+geom_vline(xintercept = lim2)
    }
    plot_par_est + theme_bw()
  })
  
  output$correlations <- renderPlot({
    df = DF()
    df = df[df$iteration %in% input$obs1:input$obs2,]
    data = df[df$rep == input$replic,2:5]
    #ggcorr(data, palette = "RdYlGn", name = "rho", 
    #       label = FALSE, label_color = "black")
    gp = ggpairs(data)
    gp
  })
  
  
  
  
  
  
  
  output$diff_fhat_hist <- renderPlot({
    if(nrow(rv$MCEM)>12){ 
      breaks <- pretty(range(rv$MCEM$diff_fhat), n = nclass.FD(rv$MCEM$diff_fhat), min.n = 1)
      bwidth <- breaks[2]-breaks[1]
      rv$MCEM$diff_fhat = c(0,diff(rv$MCEM$fhat))
      gl = ggplot(rv$MCEM[input$charts:nrow(rv$MCEM),]) + geom_histogram(aes(diff_fhat),binwidth = bwidth)#,binwidth = (max(diff_fhat)-max(diff_fhat))/50)# + ggtitle(label=paste("Last estimation:  ",mean(rv$mcem_it$K[input$charts:length(rv$mcem_it$K)])),subtitle =   paste("number of last iterations to consider: ", length(rv$mcem_it$K)-input$charts)) 
      gl  + theme_bw() + ggtitle(label = "Delta likelihood")
      #stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1))
    }
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)
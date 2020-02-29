library(emphasis)
library(GGally)
library(xtable)


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
      
      selectInput(inputId = "typePlot",
                  label = "Type of plot:",
                  choices = c("Path","Points")),
      downloadButton("report", "Generate report"),
      
      uiOutput('columns1'),
      uiOutput('columns2'),
      sliderInput("obs","Choose considered area for iterations",min=0,max=3000,value=c(0,1000)),
      
      uiOutput('reps'),
    ),
    
    
  
    
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      
     
      tableOutput("view"),
      plotOutput("parameter_estimation_general", click = "plot_click"),
      textInput("text1", "Describe the data here", value=""),
     # uiOutput("correlations_tab"),
   #   plotOutput("correlations"),
      # Output: Verbatim text for data summary ----
      #textOutput("caption"),
      
      verbatimTextOutput("summary"),
      verbatimTextOutput("code")
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  DF <- reactive({
    inFile <- input$file1
    DF_temp=NULL
    if(is.null(inFile)){
      DF_temp =     data.frame(iteration=1,par1=1,par2=1,par3=1,par4=1,E_time=1,M_time=1,sample_size=1,fhat=1,rep=1)

    }else{
     for(i in 1:nrow(inFile)){
      
        file = inFile[[i, 'datapath']]
        load(file)
      
      
        DF1 = data.frame(iteration=1:nrow(mcem),par1=mcem$par1,par2=mcem$par2,par3=mcem$par3,par4=mcem$par4,E_time=mcem$E_time,M_time=mcem$M_time,sample_size=mcem$sample_size,fhat=mcem$fhat,rep=as.character(i))
      
        DF_temp = rbind(DF_temp,DF1)
      
     }
    }
    #save(DF_temp,file="palPaper.RData")
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
  
  output$reps = renderUI({
    selectizeInput(
      inputId = "i_rep",
      label = "Select replicants:",
      choices = sort(unique(DF()$rep)),
      selected = sort(unique(DF()$rep)),
      multiple = TRUE,
      options = list(placeholder = "Begin typing title type...")
    )
  })
  
  output$columns2 = renderUI({
    selectInput(inputId = "par2",
                label = "Choose y plot parameter:",
                choices=rev(names(DF())))
  })
  
  #output$correlations_tab = renderUI({
  #  selectInput(inputId = "replic",
  #              label = "Choose replicant:",
  #              choices = unique(DF()$rep))
  #})
  

  
  output$summary <- renderPrint({
    df = DF()
    df = df[df$iteration %in% input$obs[1]:input$obs[2] & df$rep %in% input$i_rep,]
    summary(df)
  })
  

  
  output$code <- renderPrint({
    ta=NULL
    df = DF()
    df = df[df$rep %in% input$i_rep,]
    for(i in unique(df$rep)){
      ta = rbind(ta,data_to_table(df=df,replicant=i,left=input$obs[1],right=input$obs[2])) 
    }
    xtable(ta,digits=3)
  })
  

  
  output$view <- renderTable({
    ta=NULL
    df = DF()
    df = df[df$rep %in% input$i_rep,]
    for(i in unique(df$rep)){
      ta = rbind(ta,data_to_table(df=df,replicant=i,left = input$obs[1],right = input$obs[2])) 
    }
    ta
  },
  digits = 3)

  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      ta=NULL
      df = DF()
      df = df[df$rep %in% input$i_rep,]
      for(i in unique(df$rep)){
        ta = rbind(ta,data_to_table(df=df,replicant=i,left = input$obs[1],right = input$obs[2])) 
      }
      params <- list(ta = ta,text1=input$text1)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  output$parameter_estimation_general <- renderPlot({
    df = DF()
    df = df[df$rep %in% input$i_rep,]
    if(input$obs[2]>nrow(df)){
      lim2 = nrow(df)
    }else{
      lim2 = input$obs[2]
    }
    #df = df[df$iteration %in% input$obs1:lim2 & !is.na(df$mu0),]
    if(input$typePlot=="Path"){
      plot_par_est = ggplot(df,aes(colour=rep))+geom_path(aes_string(x=input$par1,y=input$par2))#+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$typePlot=="Points"){
      plot_par_est = ggplot(df,aes(colour=rep))+geom_point(aes_string(x=input$par1,y=input$par2,alpha=(1:nrow(df))/(2*nrow(df))))#+ggtitle("lambda Estimation",subtitle = paste("Current estimation:",par_est))+geom_hline(yintercept = par_est,colour="red")
    }
    if(input$par1=="iteration"){
      plot_par_est = plot_par_est +geom_vline(xintercept = input$obs[1])+geom_vline(xintercept = lim2)
    }
    plot_par_est + theme_bw()
  })
  
 # output$correlations <- renderPlot({
#    df = DF()
#   df = df[df$iteration %in% input$obs[1]:input$obs[2],]
#    data = df[df$rep == input$replic,2:5]

#    #       label = FALSE, label_color = "black")
 #   gp = ggpairs(data)
  #  gp
 # })
  
  
}

# Create Shiny app ----
shinyApp(ui, server)
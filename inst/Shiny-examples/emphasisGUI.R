

theme_emphasis =  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))  #

MCEM = data.frame(par1=NULL,par2=NULL,par3=NULL,fhat=NULL,fhat.se=NULL,E_time=NULL,M_time=NULL,mc.samplesize=NULL,effective.size=NULL,hessian.inv=NULL)

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
                             selectInput("brts_data", "Choose Phylo/Branching times:",
                                         list(
                                           "Dendroica" = "dendroica",
                                           "Anolis" = "anolis",
                                           "Cetacea" = "cetacea",
                                           "Heliconius" = "heliconius",
                                           "Plethodon" = "plethodon",
                                           "Vangidae" = "vangidae",
                                           "Foraminifera" = "foraminifera"
                                         )       
                             ),
                             h3("Initial parameters"),
                             numericInput("par1", "Initial lambda_0:", 2),
                             numericInput("par2", "Initial lambda_1:", 0.03),
                             numericInput("par3", "Initial mu_0:", 0.05),
                             
                             h3("Settings"),
                             selectInput("method", "Choose Data augmentation sampler:",
                                         list("emphasis" = "emphasis",
                                              "uniform" = "uniform"
                                         )),
                             selectInput("model", "Choose diversification model:",
                                         list("Diversity dependence" = "dd",
                                              "Exponential diversity dependence" = "edd",
                                              "Phylodiversity dependence" = "pd",
                                              "Exponential Phylodiversity dependence" = "epd"
                                         )),                
                             numericInput("sample_size", "Monte-Carlo sample size:", 10),
                             numericInput("proportion_of_subset", "Proportion of best trees to take:", 1),
                             numericInput("maxspec", "Maximum number of missing species:", 30),
                             h3("Options"),
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
                           h3("Controls"),
                           "To initialize and stop the MCEM routine.",
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
                          
                  ), #close tab panel
                  tabPanel("Analysis",
                           
                           fluidRow(
                             dataTableOutput("iterations_data"),
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
                                    textOutput("txtOutput2"),
                                    textOutput("txtOutput3")
                             )
                           )
                           
                           
                           )
                  )
                ) #close mainpanel
                
  ) # close sidebarLayout
) #close UI






server <- shinyServer(function(input,output,session) {

  input_values <- reactiveValues(brts=NULL,init_pars=NULL,brts_d=NULL)
  autoInvalidate <- reactiveTimer(intervalMs=1000,session)
  
  observe({
    # load branching times
    assign("brts",get(paste("brts_",input$brts_data,sep="")))
    input_values$brts = brts
    input_values$init_pars = c(input$par1,input$par2,input$par3)
    pars = input_values$init_pars
    autoInvalidate() 
  })
  
  observeEvent(input$gogobutt, { isolate({ rv$run=T      }) })
  observeEvent(input$stopbutt, { isolate({ rv$run=F      }) })
  
  output$phylogenetic_tree <- renderPlot({
    tree = list(wt = -diff(c(input_values$brts,0)),to=rep(1,length(input_values$brts)-1))
    tree = vectors2phylo(tree)
    plot(tree,show.tip.label=F,root.edge=TRUE)
    axisPhylo()
  })
  

}
)

shinyApp(ui, server)
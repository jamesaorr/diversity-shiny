#######################################
############# Libraries ###############
#######################################

library(shiny)
library(vegan)
library(ggplot2)
library(bslib)

#######################################
########### User interface ############
#######################################

ui <- fluidPage(
  withMathJax(),
  theme = bs_theme(bootswatch = "flatly"),
  titlePanel("Diversity Metrics"),
  
  tags$style(HTML("
    body {
      margin: 40px;
    }
    .justified-text {
      text-align: justify;
    }
  ")),
  
  fluidRow(
    column(6,
           tags$p(class = "justified-text",
                  "This Shiny app lets you explore the most common diversity metrics that are used in ecology. 
                  You can enter species abundances in the text box below and the app will produce the corresponding rank abundance plot. 
                  The app will also print the values of the traditional diversity metrics for your community. 
                  Additionally, you can explore Hill diversity, which connects all of the traditional metrics and gives them all the same units."
           ),
           br(),
           textInput("species_abundances", "Enter species abundances:", 
                     "50, 2, 5, 6, 4, 8, 10, 1, 6, 1")
    ),
    column(5, plotOutput("abundance_curve"))
  ),
  
  br(),
  
  tabsetPanel(
    br(), 
    tabPanel("Traditional Diversity Metrics",
             tags$h3("Traditional Diversity Metrics"),
             br(),
             fluidRow(
               column(3, wellPanel(
                 style = "display: flex; flex-direction: column; justify-content: center; align-items: center; height: 310px;",
                 tags$h4("Richness"),
                 tags$p("$$S$$"),
                 tags$p("The number of species in a community"),
                 verbatimTextOutput("richness_out2")
               )),
               
               column(3, wellPanel(
                 style = "display: flex; flex-direction: column; justify-content: center; align-items: center; height: 310px;",
                 tags$h4("Shannon index"),
                 tags$p("$$H^{\\prime}=-\\sum_{i=1}^S p_i \\ln \\left(p_i\\right)$$"),
                 tags$p("The uncertainty in predicting the identity of a randomly chosen individual"),
                 verbatimTextOutput("shannon_out")
               )),
               
               column(3, wellPanel(
                 style = "display: flex; flex-direction: column; justify-content: center; align-items: center; height: 310px;",
                 tags$h4("Shannon evenness"),
                 tags$p("$$J^{\\prime} = \\frac{H^{\\prime}}{\\ln(S)}$$"),
                 tags$p("A unitless index of how evenly individuals are distributed among species"),
                 verbatimTextOutput("evenness_out")
               )),
               
               column(3, wellPanel(
                 style = "display: flex; flex-direction: column; justify-content: center; align-items: center; height: 310px;",
                 tags$h4("Simpson Index"),
                 tags$p("$$D=\\sum_{i=1}^R p_i^2$$"),
                 tags$p("The probability that two individuals will belong to the same species"),
                 verbatimTextOutput("simpson_out")
               ))
             )
    ),
    
    br(), 
    
    tabPanel("Hill Diversity Metrics",
             tags$h3("Hill Diversity Metrics"),
             fluidRow(
               column(4, 
                      wellPanel(
                        tags$h4("Effective Diversity & Hill Diversity"),
                        tags$p("Effective diversity refers to the number of equally abundant species that would result in the same value of that diversity metric (e.g., Shannon or Simpson). 
                               It is helpful because the units, rather than being in information or probability, become the number of species. 
                               Weirdly, we can now have non-integer numbers of species!"),
                        tags$p("Hill diversity is a family of metrics that generalizes all of the traditional metrics into one framework. 
                                Through Hill diversity we can see that metrics like richness, Shannon, or Simpson put different weight on rare vs. common species.
                                The parameter 'q' controls the sensitivity to rare/common species, with q=0 corresponding to richness, q=1 to Shannon, and q=2 to Simpson.")
                      )
               ),
               
               column(8, 
                      
                      fluidRow(
                        column(4, wellPanel(
                          tags$h4("Richness"),
                          tags$p("$$S$$"),
                          verbatimTextOutput("richness_out")
                        )),
                        
                        column(4, wellPanel(
                          tags$h4("Effective Shannon"),
                          tags$p("$$exp(H^{\\prime})$$"),
                          verbatimTextOutput("effective_shannon_out")
                        )),
                        
                        column(4, wellPanel(
                          tags$h4("Effective Simpson"),
                          tags$p("$$\\frac{1}{D}$$"),
                          verbatimTextOutput("inv_simpson_out")
                        ))
                      ),
                      
                      br(),
                      
                      fluidRow(
                        column(6, 
                               wellPanel(
                                 tags$h4("Hill Diversity"),
                                 tags$p("$$D_q = (\\sum p_i^q)^{1/(1-q)}$$"),
                                 verbatimTextOutput("hill_out")
                               )
                        ),
                        column(6, 
                               wellPanel(
                                 sliderInput("hill_q", "Hill exponent (q):", min = 0, max = 2, value = 1, step = 0.1)
                               )
                        )
                      ),
               )
             )
    )
  )
)


#######################################
############## Server #################
#######################################

server <- function(input, output) {
  parsed_abundances <- reactive({
    # Remove trailing commas and extra spaces before splitting
    cleaned_input <- gsub(",\\s*$", "", input$species_abundances)
    as.numeric(unlist(strsplit(cleaned_input, ",")))
  })
  
  output$abundance_curve <- renderPlot({
    abundances <- parsed_abundances()
    abundances <- abundances[abundances > 0]  # Remove zeros
    df <- data.frame(Species = seq_along(abundances), Abundance = sort(abundances, decreasing = TRUE))
    
    ggplot(df, aes(x = Species, y = Abundance)) +
      geom_bar(stat = "identity", fill = "darkgoldenrod1") +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 14),  
        axis.title = element_text(size = 16),  
        axis.title.x = element_text(margin = margin(t = 10)),  
        axis.title.y = element_text(margin = margin(r = 10))   
      ) +
      labs(x = "Rank", y = "Abundance")
  })
  
  output$richness_out <- renderText({
    abundances <- parsed_abundances()
    abundances <- abundances[abundances > 0]  # Remove zeros
    length(abundances)
  })
  
  # have to put it in twice as I call it twice (on each tab)
  output$richness_out2 <- renderText({
    abundances <- parsed_abundances()
    abundances <- abundances[abundances > 0]  # Remove zeros
    length(abundances)
  })
  
  output$shannon_out <- renderText({
    round(diversity(parsed_abundances(), index = "shannon"), 3)
  })
  
  output$evenness_out <- renderText({
    abundances <- parsed_abundances()
    abundances <- abundances[abundances > 0]  # Remove zeros
    S <- length(abundances)
    shannon <- diversity(abundances, index = "shannon")
    round(shannon / log(S), 3)
  })
  
  output$simpson_out <- renderText({
    round(1 - diversity(parsed_abundances(), index = "simpson"), 3)
  })
  
  output$effective_shannon_out <- renderText({
    round(exp(diversity(parsed_abundances(), index = "shannon")), 3)
  })
  
  output$inv_simpson_out <- renderText({
    round(diversity(parsed_abundances(), index = "invsimpson"), 3)
  })
  
  output$hill_out <- renderText({
    round(renyi(parsed_abundances(), hill = TRUE, scales = input$hill_q), 3)
  })
}

shinyApp(ui, server)

library(shiny)

renderInputs <- function(prefix) {
  wellPanel(
    fluidRow(
      column(6,
             h5("Your Data"),     
             numericInput(paste0(prefix, "_", "samplesize"), "Sample Size:", 200, min=1),
             
             checkboxInput(paste0(prefix, "_", "customize_data"), "Additional Custom Features",value=FALSE),
             conditionalPanel(
               condition =paste0("input.",prefix, "_", "customize_data"),
               selectInput(paste0(prefix, "_", "missing"), "Missing Data:",
                           c("None" = FALSE,
                             "At Random"=TRUE))),
             h5("Methods"),
             selectInput(paste0(prefix, "_", "rotation"), "Rotation:",
                         c("Oblimin"="oblimin","Varimax"="varimax", "Promax"="promax","None"="none")),
             selectInput(paste0(prefix, "_", "method"), "Factoring Method:",
                         c("Maximum Likelihood FA"="mle",
                           "Principal Axis Factor Analysis"="pa",
                           "minimum residual (OLS) factoring"="minres" ,
                           "Principal Components"="pc" ))
      ),
      column(6,
             h5("True Model"),
             sliderInput(paste0(prefix, "_", "nfactors"), "Number of Factors:", 
                         min=0, max=10, value=5),
             sliderInput(paste0(prefix, "_", "rfactors"), "Correlation between Factors:", 
                         min=-1, max=1, value=.05,step=.05),
             
             sliderInput(paste0(prefix, "_", "loading"), "Factor Loading:", 
                         min=0, max=1, value=.5,step=.05),
             numericInput(paste0(prefix, "_", "items_p_f"), "Number of items per Factor", min=1, value=5),
             numericInput(paste0(prefix, "_", "itemsR_p_f"), "How many items are reverse scored per Factor?", min=0, value=0),
             checkboxInput(paste0(prefix, "_", "customize_model"), "Additional Custom Features",value=FALSE),
             
             conditionalPanel(condition = paste0("input.",prefix, "_", "customize_model"),
                              checkboxInput(paste0(prefix, "_", "loading_norm"), "Customize Factor Loadings",value=FALSE),
                              conditionalPanel(
                                condition = paste0("input.",prefix, "_", "loading_norm"),
                                numericInput(paste0(prefix, "_", "loading_norm_sd"), "Standard Deviation for Factor Loadings", min=0, value=0,max=1)),
                              checkboxInput(paste0(prefix, "_", "r_norm"), "Customize Correlations between Factors",value=FALSE),
                              conditionalPanel(
                                condition = paste0("input.",prefix, "_", "r_norm"),
                                numericInput(paste0(prefix, "_", "r_norm_sd"), "Standard Deviation for Correlations ", min=0, value=0,max=1)),
                              checkboxInput(paste0(prefix, "_", "custom_item"), "Customize distribution of items by factor",value=FALSE),
                              conditionalPanel(
                                condition = paste0("input.",prefix, "_", "custom_item"),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 0 "),
                                  numericInput(paste0(prefix, "_", "f1_items"), "Number of items for Factor 1", min=1, value=5),
                                  numericInput(paste0(prefix, "_", "f1_itemsR"), "How many Factor 1 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 1 "),
                                  numericInput(paste0(prefix, "_", "f2_items"), "Number of items for Factor 2", min=1, value=5),
                                  numericInput(paste0(prefix, "_", "f2_itemsR"), "How many Factor 2 items are reverse scored?", min=0,value=0,max="input.f2_items")),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 2 "),
                                  numericInput(paste0(prefix, "_", "f3_items"), "Number of items for Factor 3", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f3_itemsR"), "How many Factor 3 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 3 "),
                                  numericInput(paste0(prefix, "_", "f4_items"), "Number of items for Factor 4", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f4_itemsR"), "How many Factor 4 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 4 "),
                                  numericInput(paste0(prefix, "_", "f5_items"), "Number of items for Factor 5", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f5_itemsR"), "How many Factor 5 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 5 "),
                                  numericInput(paste0(prefix, "_", "f6_items"), "Number of items for Factor 6", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f6_itemsR"), "How many Factor 6 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 6 "),
                                  numericInput(paste0(prefix, "_", "f7_items"), "Number of items for Factor 7", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f7_itemsR"), "How many Factor 7 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 7 "),
                                  numericInput(paste0(prefix, "_", "f8_items"), "Number of items for Factor 8", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f8_itemsR"), "How many Factor 8 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 8 "),
                                  numericInput(paste0(prefix, "_", "f9_items"), "Number of items for Factor 9", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f9_itemsR"), "How many Factor 9 items are reverse scored?", min=0, value=0)),
                                conditionalPanel(
                                  condition = paste0("input.",prefix, "_", "nfactors > 9 "),
                                  numericInput(paste0(prefix, "_", "f10_items"), "Number of items for Factor 10", min=1,  value=5),
                                  numericInput(paste0(prefix, "_", "f10_itemsR"), "How many Factor 10 items are reverse scored?", min=0, value=0))))
             ,
             h5("Simulation"),
             numericInput(paste0(prefix, "_", "seed"), "Seed:", 12345),
             sliderInput(paste0(prefix, "_", "ndatasets"), "Number of replications:", min = 1, max = 1000, value = 200)
             
      )
    ),
    p(actionButton(paste0(prefix, "_", "recalc"),
                   "Re-run simulation", icon("random")
    ))
  )
}

# Define UI for application that plots random distributions
fluidPage(theme="simplex.min.css",
          tags$style(type="text/css",
                     "label {font-size: 12px;}",
                     ".recalculating {opacity: 1.0;}"
          ),
          
          # Application title
          tags$h2("Factor Enumeration: Creating Custom Model Fit Criteria with Simulations"),
          p("By S. Mason Garrison"),p("More information available",
                                      tags$a(href="https://www.dropbox.com/s/52s1t27jw4wzs6d/arp2017.pdf?dl=0", "here (ARP Poster)"),
                                      "and",
                                      tags$a(href="https://www.dropbox.com/s/vpflvbf4vd3s2nc/talk.pdf?dl=0", "here (Talk)"),
                                      "."),
          hr(),
          
          fluidRow(
            column(6, 
                   tags$h3("Simulation A")),
            column(6, tags$h3("Simulation B"))
          ),
          fluidRow(
            column(6, 
                   renderInputs("a")),
            column(6, renderInputs("b"))
          ),
          fluidRow(
            column(6,
                   print(resultsA)
            ),
            column(6,
                   print(resultsB)
            )
          )
)

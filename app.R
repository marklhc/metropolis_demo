#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
# par(mar = c(4, 4, 2, 0) + 0.1)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Metropolis Algorithm"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("sd", "Proposal SD", 0.05, 0.5, value = 0.3),
      actionButton("init", "Initialize"), 
      actionButton("prop", "Propose a value"), 
      actionButton("acc", "Accept/Reject"), 
      actionButton("run100", "Draw 100 samples"), 
      textOutput("nsample"), 
      textOutput("acc_rate")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", 
                 plotOutput("TargetPlot"), 
                 plotOutput("SamplePlot")), 
        tabPanel("Trace", 
                 plotOutput("TracePlot"), 
                 plotOutput("AcfPlot"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  dens_kern <- function(x) {
    dnorm(x, 0.5, 1 / sqrt(2)) * dbeta(x, 13, 9)
  }
  x0 <- seq(-0.3, 0.3, length.out = 101)
  v <- reactiveValues(th0 = NULL, proposed = NULL, sam = NULL, 
                      accept = 0)
  observeEvent(input$init, {
    v$th0 <- runif(1)
    v$proposed <- NULL
    v$sam <- NULL
    v$accept <- 0
  })
  observeEvent(input$prop, {
    if (!is.null(v$th0)) {
      v$proposed <- rnorm(1, sd = input$sd) + v$th0
    }
  })
  observeEvent(input$acc, {
    if (!is.null(v$proposed)) {
      u <- runif(1)
      if (u < dens_kern(v$proposed) / dens_kern(v$th0)) {
        v$th0 <- v$proposed
        v$accept <- v$accept + 1
      }
      v$sam <- c(v$sam, v$th0)
      v$proposed <- NULL
    }
  })
  observeEvent(input$run100, {
    if (!is.null(v$th0)) {
      proposed_jump <- rnorm(100, sd = input$sd)
      sam_mcmc <- rep(NA, length = 101)
      sam_mcmc[1] <- v$th0
      for (i in 1:100) {
        proposed <- sam_mcmc[i] + proposed_jump[i]
        if (runif(1) < dens_kern(proposed) / dens_kern(sam_mcmc[i])) {
          sam_mcmc[i + 1] <- proposed
          v$accept <- v$accept + 1
        } else {
          sam_mcmc[i + 1] <- sam_mcmc[i]
        }
      }
      v$th0 <- sam_mcmc[101]
      v$sam <- c(v$sam, sam_mcmc[2:101])
      v$proposed <- NULL
    }
  })
  output$TargetPlot <- renderPlot({
    if (!is.null(v$th0)) {
      prop_sd <- input$sd
      curve(dens_kern(x), from = 0, to = 1, 
            xlab = expression(theta), ylab = "")
      curve(dnorm(x, v$th0, prop_sd) * 0.1, from = v$th0 - 0.3, 
            to = v$th0 + 0.3, add = TRUE, col = "lightblue")
      x <- x0 + v$th0
      polygon(c(v$th0 - 0.3, x, v$th0 + 0.3), 
              c(0, dnorm(x, v$th0, prop_sd) * 0.1, 0), 
              col = rgb(0.8, 1, 1, 0.2), border = NA)
      points(v$th0, 0, pch = 19, col = "skyblue3")
      segments(v$th0, 0, y1 = dens_kern(v$th0), col = "blue")
      if (!is.null(v$proposed)) {
        proposed <- v$proposed
        segments(proposed, 0, y1 = dens_kern(proposed), col = "red")
        points(proposed, 0, pch = 19, col = "red")
        text(0, 2, bquote(italic(P)(Accept) == 
                            .(round(min(dens_kern(v$proposed) / 
                                          dens_kern(v$th0), 1), 3) * 100) ~ "%"), 
             pos = 4)
      }
    }
  })
  output$SamplePlot <- renderPlot({
    if (!is.null(v$sam)) {
      hist(v$sam, col = "skyblue3", xlab = expression(theta), 
           xlim = c(0, 1), main = "Sampled Values", freq = FALSE)
      if (length(v$sam) > 10) {
        lines(density(v$sam, bw = "SJ"), col = "blue")
      }
    }
  })
  output$TracePlot <- renderPlot({
    if (!is.null(v$sam)) {
      plot(v$sam, type = "l", xlab = "Sample", ylab = "")
    }
  })
  output$AcfPlot <- renderPlot({
    if (!is.null(v$sam)) {
      acf(v$sam, main = "Autocorrelation Plot")
    }
  })
  output$nsample <- renderText({ 
    paste("Number of simulation samples:", length(v$sam))
  })
  output$acc_rate <- renderText({ 
    paste("Acceptance rate:", round(v$accept / length(v$sam), 2))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)


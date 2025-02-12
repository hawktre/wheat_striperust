library(shiny)
library(circular)
library(BH)
library(tidyverse)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Kernel Exploration"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("mu",
                  "mu (mean of f)",
                  min = 0,
                  max = 2*pi,
                  value = pi),
      sliderInput("nu",
                  "nu (mean of g)",
                  min = 0,
                  max = 2*pi,
                  value = pi),
      sliderInput("delta",
                  "delta (dispersion of f)",
                  min = 0.01,
                  max = 20,
                  value = 1.25),
      sliderInput("kappa",
                  "kappa (dispersion of g)",
                  min = 0.01,
                  max = 20,
                  value = 3),
      sliderInput("g0",
                  "g0 (distance scale)",
                  min = 1,
                  max = 10,
                  value = 3),
      sliderInput("b",
                  "b (power law exponent)",
                  min = 2.01,
                  max = 3.00,
                  value = 2.25)
    ),
    
    
    mainPanel(
      tabsetPanel(
        tabPanel("Kernel Contours", plotOutput('kernPlot')),
        tabPanel("f, g Function Plot", plotOutput('fgPlot')),
        tabPanel("Histogram of Raw Kernel Values", plotOutput('scalePlot'))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  f_fn <- function(x){
    dvonmises(circular(x), 
              circular(input$mu), 
              (input$delta))
  }
  
  g_fn <- function(x){
    dvonmises(circular(x), 
              circular(input$nu), 
              (input$kappa))*input$g0
  }
  
  K_fn <- function(rho, phi){
    f_fn(phi)*(input$b - 1)*(input$b - 2)/(g_fn(phi)^2)*((1 + rho/g_fn(phi))^(-input$b))
  }
  
  rho_vals <- seq(1, 10, length = 50)
  phi_vals <- seq(0, 2*pi, length = 50)
  
  output$fgPlot <- renderPlot({
    ylim_upr <- max(max(g_fn(seq(from = 0, to = 2*pi, length = 100))/input$g0),
                    max(f_fn(seq(from = 0, to = 2*pi, length = 100)))
    )
    
    curve(f_fn, from = 0, to = 2*pi, col = 1,
          ylab = expression(paste('fn(', theta, ')', sep = '')),
          xlab = expression(theta),
          ylim = c(0, ylim_upr))
    curve(g_fn(x)/input$g0, from = 0, to = 2*pi, add = T, col = 2)
    legend(x = 'topright', legend = c('f', 'g/g0'), col = 1:2, lty = 1)
  })
  
  output$kernPlot <- renderPlot({
    
    grid_df <- mutate(expand.grid(rho_vals, phi_vals),
                      z = K_fn(Var1, Var2))
    
    
    ggplot(grid_df, aes(x = Var1, y = Var2, z = (log(z)))) + 
      geom_contour(binwidth = 0.5, 
                   aes(color = stat(level))) +
      scale_y_continuous(breaks = seq(0, 2*pi, length = 5), 
                         limits = c(0, 2*pi),
                         labels = expression(0, pi/2, pi, 3*pi/2, 2*pi)) +
      xlim(range(rho_vals)) +
      scale_color_continuous(limits = c(-10, -1)) +
      labs(x = '', y = '') +
      guides(color = guide_colorbar('log(K(x, y))')) +
      coord_polar(theta = "y") +
      theme_bw()
    
  })
  
  output$scalePlot <- renderPlot({
    grid_df <- mutate(expand.grid(rho_vals, phi_vals),
                      z = K_fn(Var1, Var2))
    
    ggplot(grid_df, aes(x = z)) +
      geom_histogram(bins = 35) +
      theme_bw() +
      labs(x = 'K(x, y)', y = 'count')
  })
}

# Run the application 
shinyApp(ui = ui, server = server)


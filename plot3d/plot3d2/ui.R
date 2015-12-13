options(rgl.useNULL=TRUE)

library(shinyRGL)
library(rgl)

shinyUI(fluidPage(title="K29 PCA Training Dataset", 
          sidebarLayout(
            sidebarPanel(
              h4("Training Data Controls"),
              sliderInput("tpts", 
                          "Number of background points:", 
                          min = 1000, 
                          max = 5686, 
                          value = 1000),
              hr(),    
              sliderInput("tbins", 
                          "Number of bins:", 
                          min = 1, 
                          max = 73, 
                          step = 1,
                          value = 1),
              hr(),    
              checkboxInput("textflag", "Bin Labels", FALSE),
              helpText(HTML("Created using <a href = \"http://github.com/trestletech/shinyRGL\">shinyRGL</a>."))),
            mainPanel(
              webGLOutput("myPlotRGLT", width="800px", height="800px")
            )
        )
))
        

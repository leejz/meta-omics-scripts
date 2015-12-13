options(rgl.useNULL=TRUE)

library(shinyRGL)
library(rgl)

shinyUI(fluidPage(title="k29 PCA Full Dataset",
          sidebarLayout(
            sidebarPanel(
              h4("Full Dataset Controls"),
              hr(),
              sliderInput("fpts", 
                          "Number of background points:", 
                          min = 1000, 
                          max = 51562, 
                          value = 1000),
              hr(),    
              sliderInput("fbins", 
                          "Number of bins:", 
                          min = 1, 
                          max = 73, 
                          step = 1,
                          value = 1),
              hr(),    
              checkboxInput("textflag", "Bin Labels", FALSE),
              helpText(HTML("Created using <a href = \"http://github.com/trestletech/shinyRGL\">shinyRGL</a>."))),
            mainPanel(
              webGLOutput("myPlotRGLF", width="800px", height="800px")
            )
        )
))
        

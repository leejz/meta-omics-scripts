# Must be executed BEFORE rgl is loaded on headless devices.
options(rgl.useNULL=TRUE)

library(shiny)
library(shinyRGL)
library(rgl)
library(RColorBrewer)

d3 <- read.csv("d3.k29.csv", sep = " ")

shinyServer(function(input, output) {
  #zoom <- 1
  #userMatrix <- matrix(c(1, 0, 0, 0, 0, 0.3420201, 0.9396926, 0, 0, -0.9396926, 0.3420201, 0, 0, 0, 0, 1), ncol=4,nrow=4)
  #indowRect <-c(0, 0, 256, 256)
  
  output$myPlotRGLF <- renderWebGL({
    bintext <-{}
    d3_unique_cluster <- unique(d3$svm.pred[d3$svm.pred != 0])
    bkxPts <- d3$HPminus
    bkyPts <- d3$HPplus
    bkzPts <- d3$HP3
    
    d3subset <- d3[d3$svm.pred != 0 & d3$svm.pred %in% d3_unique_cluster[1:input$fbins],]
    fgxPts <- d3subset$HPminus
    fgyPts <- d3subset$HPplus
    fgzPts <- d3subset$HP3
    fgcol <- d3_unique_cluster[d3subset$svm.pred]
      
    # Add a label in the center of each extracted bin
    dbspan <- length(d3_unique_cluster)
    gbr<-colorRampPalette(c("green","blue","orange","red"))        
    palette(adjustcolor(gbr(dbspan+1), alpha.f = 1))
    for (i in 1:length(d3_unique_cluster)) {
      bintext <- rbind(bintext, c(as.character(d3_unique_cluster[i]),
                                  mean(d3$HPminus[d3$svm.pred == d3_unique_cluster[i]]), 
                                  mean(d3$HPplus[d3$svm.pred == d3_unique_cluster[i]]), 
                                  mean(d3$HP3[d3$svm.pred == d3_unique_cluster[i]])))
    }
 
    #view in 3d  
    #open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
    #par3d(cex=.6)
    plot3d(bkxPts[1:input$fpts],
         bkyPts[1:input$fpts],
         bkzPts[1:input$fpts],
         col=rgb(0,0,0), 
         size=3, 
         type='p',
         alpha=0.3,
         xlab="PC1",
         ylab="PC2",
         zlab="PC3")
  
    points3d(fgxPts,
         fgyPts,
         fgzPts,
         col=fgcol+1,
         size=3,
         alpha=.75)
  
    if (input$textflag == TRUE) {
      text3d(x=bintext[,2], y=bintext[,3], z=bintext[,4],text = bintext[,1])
    }
  
    axes3d()
    zoom<-par3d()$zoom
    userMatrix<-par3d()$userMatrix
    windowRect<-par3d()$windowRect  
  })
  
})
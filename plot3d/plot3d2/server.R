# Must be executed BEFORE rgl is loaded on headless devices.
options(rgl.useNULL=TRUE)

library(shiny)
library(shinyRGL)
library(rgl)
library(RColorBrewer)

d2 <- read.csv("d2.k29.csv", sep = " ")

shinyServer(function(input, output) {
  #zoom <- 1
  #userMatrix <- matrix(c(1, 0, 0, 0, 0, 0.3420201, 0.9396926, 0, 0, -0.9396926, 0.3420201, 0, 0, 0, 0, 1), ncol=4,nrow=4)
  #windowRect <-c(0, 0, 256, 256)
  
  output$myPlotRGLT <- renderWebGL({
    bintext <-{}
    d2_unique_cluster <- unique(d2$db_cluster[d2$db_cluster != 0])
    bkxPts <- d2$HPminus
    bkyPts <- d2$HPplus
    bkzPts <- d2$HP3
    
    d2subset <- d2[d2$db_cluster != 0 & d2$db_cluster %in% d2_unique_cluster[1:input$tbins],]
    fgxPts <- d2subset$HPminus
    fgyPts <- d2subset$HPplus
    fgzPts <- d2subset$HP3
    fgcol <- d2_unique_cluster[d2subset$db_cluster]
    
    # Add a label in the center of each extracted bin
    dbspan <- length(d2_unique_cluster)
    gbr<-colorRampPalette(c("green","blue","orange","red"))    
    palette(adjustcolor(gbr(dbspan+1), alpha.f = 1))
    for (i in 1:length(d2_unique_cluster)) {
      bintext <- rbind(bintext, c(as.character(d2_unique_cluster[i]),
                                  mean(d2$HPminus[d2$db_cluster == d2_unique_cluster[i]]), 
                                  mean(d2$HPplus[d2$db_cluster == d2_unique_cluster[i]]), 
                                  mean(d2$HP3[d2$db_cluster == d2_unique_cluster[i]])))
    }
  
  #view in 3d  
  #open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
  #par3d(cex=.6)
  plot3d(bkxPts[1:input$tpts],
         bkyPts[1:input$tpts],
         bkzPts[1:input$tpts],
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
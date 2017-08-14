
SPIplots = function(MergeFile){
  A = arrange(MergeFile, desc(Best.Z)) #arrange largest to smallest
  lim = ceiling(max(A$Best.Z)) # find best limits for graph
  plot(A$Best.Z, xlab = "GFP strains", ylab = " Z-score", pch = 19,col = rainbow(8000),cex = .4,panel.first=grid(20,10),ylim = c(-1*lim,lim))
  abline(h = c(2,-2)) # add lines at ±2
  
  
  B = filter(MergeFile, (GOI.Z != "n/a")&(GBP.Z != "n/a"))
  Xlim = ceiling(as.numeric(max(B$GOI.Z))) # Best limits for graph
  Ylim = ceiling(as.numeric(max(B$GBP.Z)))
  Color = 9*(B$Best.Z>=2)+1 # Vector of colors, red = SPI, black  = non-SPI
  plot(B$GOI.Z,B$GBP.Z,pch = 19,panel.first=grid(12,12),asp =1,cex = .4, xlim = c(-1*Xlim,Xlim),ylim = c(-1*Ylim,Ylim), xlab = "GOI control Z-score", ylab = "GBP control Z-score", col = Color, main = "Control Correlation", )
  abline(h = 2,v=2)
  lines(c(100,-96),c(-96,100)) # Line showing (x+y)/2 = 2, line showing SPI limit
  legend("bottomleft",c("SPI","non-SPI"), col = c(10,1), pch =19)
}

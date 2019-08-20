
MMfit = function(SDataFile,componentPlot = F, fitPlot = F){
  
  # This is a basic script to apply a mixture model and calculate q(x) the probability of inclusion in the hit distribution for a set of data points,
  # it can also plot the components of the distribution and the fit of the model to the data. 
  is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
  if(!is.installed("tidyverse")){
    install.packages("tidyverse")
  }
  if(!is.installed("mclust")){
    install.packages("mclust")
  }
  if(!is.installed("ggplot2")){
    install.packages("ggplot2")
  }
  library(tidyverse)
  library(mclust)
  library(ggplot2)
  
  # This script takes in a simple Name Value .csv file and creates a mixture model and calculates the q(x) values
  select = dplyr::select
  # Prepare files
  sData = read_csv(SDataFile)
  # Fit models
  modG2 = densityMclust(sData$Value, G = 2)
  
  # The order of the fit models is not always such that m_1<m_2 so define ordering based on this
  Order = sort(modG2$parameters$mean, index.return = T)$ix
  maxG = 2
  minG = 1
  means = modG2$parameters$mean[Order]
  pros = modG2$parameters$pro[Order]
  sds = sqrt(modG2$parameters$variance$sigmasq[Order])
  
  # Define whether the model fitting is successful based on whether the means are within 1.5 sd of each other
  failure = 0
  if((means[2]- means[1])<1.5*sds[1]){
    failure = 1
  }
  print(failure)
  if(failure==0){
    
  
  # Define distributions: f1 is mixture model, f2 is hit peak component
  f1 = function(x) modG2$parameters$pro[minG]*dnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]))+modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
  f2 = function(x) modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
  
  # Work out q(x) and write data
  
  q = function (x) f2(x)/f1(x)
  
  sData2 = mutate(sData, qval = q(Value))
  
  if(substr(SDataFile,start = nchar(SDataFile)-4,stop = nchar(SDataFile))==".csv"){
    name = paste0(substr(SDataFile,start = 1,stop = nchar(SDataFile)-4), "MMqval.csv")
  } else{name = paste0(SDataFile,"MMqval.csv")}
  write_csv(sData2,name)

  # Plots data
  
  x = seq(min(sData$Value),max(sData$Value),length = 1000)
  dist = f1(x)
  if(fitPlot){
    data = hist(sData$Value, breaks = 150, plot = FALSE)
    distVals = tibble(X = x, Y = dist, Distribution = "Mixture model")
    dataVals = tibble(X = data$mids, Y = data$density, Distribution = "Screen data")
    Vals = rbind(distVals, dataVals)
    png(filename = paste0(name,"MMfit.png"), width = 600)
    print(ggplot(Vals, aes(x= X, y = Y, colour = Distribution))+ geom_line(linetype = 1) +  theme_bw() +scale_color_brewer(palette="Set1")+labs(x="LGR",y="density")+theme(text = element_text(size = 24,family = "Helvetica")) 
    )
    dev.off()}
  
  if(componentPlot){
    comp = tibble(X =x, Y = modG2$parameters$pro[1]*dnorm(x, mean = modG2$parameters$mean[1], sd = sqrt(modG2$parameters$variance$sigmasq[1])), Distribution = "Distribution 1")
    int= tibble(X=x, Y = modG2$parameters$pro[2]*dnorm(x, mean = modG2$parameters$mean[2], sd = sqrt(modG2$parameters$variance$sigmasq[2])), Distribution = paste("Distribution", 2))
    comp = rbind(comp, int)
    comp = rbind(comp, tibble(X=x, Y=dist, Distribution = "Mixture model"))
    png(filename = paste0(name,"MMcomponents.png"), width = 600)
    print(ggplot(comp, aes(x = X, y = Y, colour = Distribution)) + geom_line() +  theme_bw() +scale_color_brewer(palette="Set1")+labs(x="LGR",y="density")+theme(text = element_text(size = 24,family = "Helvetica"),legend.background = element_rect(fill="transparent"))
    )
    dev.off()
    
  }

  return(list(means=means, pros=pros, sds=sds))
  } else {print("Mixture model fitting failed")}
}


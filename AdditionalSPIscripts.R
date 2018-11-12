# Version 1.0 15/10/18
# Note that this code should use files output from AnalyseSPIdata or derivatives to work

removeOutliers = function(Data){
  # based on Tukey [Tukey, John W (1977). Exploratory Data Analysis. Addison-Wesley]
  # This is used in summariseValiData
  Q1 = quantile(Data,1/4)
  Q3 = quantile(Data,3/4)
  Range = c(Q1-1.5*(Q3-Q1), Q3 + 1.5*(Q3-Q1))
  
  modData = Data[which((Data>Range[1])&(Data<Range[2]))]
  return(modData)
}


writeHits = function(cutoff, sData, screen, name){
  # This code can be used to write a file of Hits in basic text format, useful for GO analysis
  Hits = filter(sData, Best.LGR>cutoff)
  write_csv(tibble(ORF =Hits$ORF), paste0(screen, name,"Hits.csv"))
  write_csv(tibble(ORF=sData$ORF), paste0(screen,"Background.csv"))
  
}

summariseValiData = function(SDataFile, VDataFile, writeFile){
  # This code creates a file with comparison of the LGR from the screen and in the validation
  # It also calculates validation cutoffs as mean+2*sd of controls and includes a Boolean "Validated?" column
  select = dplyr::select
  # Prepare files
  sData = read_csv(SDataFile)
  vData = read_csv(vDataFile)
  ctrls = filter(vData, ORF %in% c("Positive control", "Positive Control", "positive control"))
  vDataF = filter(vData, !(ORF %in% c("Positive control", "Positive Control", "positive control")))
  
  # Prepare cutoff for each plate
  NvalPlates = range(ctrls$Plate)[2]
  Vcutoff =c()
  for(k in 1:NvalPlates){
    Vctrls = removeOutliers(filter(ctrls,Plate ==k)$Best.LGR)
    #Vctrls = filter(ctrls,Plate ==k)$Best.LGR
    Vcutoff[k] = mean(Vctrls)+2*sd(Vctrls)
  }
  
  # Arrange data
  sData1 = select(sData, ORF, sLGR = Best.LGR, sZ = Best.Z)
  vDataF1 = mutate(vDataF, Vcutoff = Vcutoff[vDataF$Plate], `Validated?` = Best.LGR>Vcutoff)
  vDataF1 = select(vDataF1, ORF, Plate, vLGR = Best.LGR, Vcutoff, `Validated?`)
  allData = inner_join(sData1, vDataF1, by = "ORF")
  write_csv(allData, writeFile)
  
}

DPmodelfitBasic = function(SDataFile,Screen,compDistPlot = F, fitPlot = F){
  # The basic version does not include a validation screen
  select = dplyr::select
  # Prepare files
  sData = read_csv(SDataFile)
  # Fit models
  modG2 = densityMclust(sData$Best.LGR, G = 2)
  
  Zcutoff = mean(sData$Best.LGR)+2*sd(sData$Best.LGR)
  
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
  # Define distributions: f1 is mixture model, f2 is hit peak component
  f1 = function(x) modG2$parameters$pro[minG]*dnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]))+modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
  f2 = function(x) modG2$parameters$pro[maxG]*dnorm(x, mean = modG2$parameters$mean[maxG], sd = sqrt(modG2$parameters$variance$sigmasq[maxG]))
  
  # Find all cutoffs
  
  q = function (x) f2(x)/f1(x)
  interceptG2q5s =c()
  interceptG2q50s =c()
  for(k in 1:10){ # Optim can be quite problematic, only use on small regions then find most reasonable fit over all regions
    #interceptG2ps[k]= (optim(c(0),f = function(x) (f2(x)-0.95*(f1(x)))^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    interceptG2q5s[k]= (optim(c(0),f = function(x) (q(x)-0.95)^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    interceptG2q50s[k]= (optim(c(0),f = function(x) (q(x)-0.5)^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    if((interceptG2q5s[k] > (((k+1)*0.1)-0.049))|(interceptG2q5s[k] < (((k-1)*0.1)+0.049))){
      interceptG2q5s[k] = 0
    }
    if((interceptG2q50s[k] > (((k+1)*0.1)-0.049))|(interceptG2q50s[k] < (((k-1)*0.1)+0.049))){
      interceptG2q50s[k] = 0
    }
  }
  Xq5 = max(interceptG2q5s)
  Xq50 = max(interceptG2q50s)

  interceptG2ps =c()
  interceptG2p1s =c()
  
  for(k in 1:10){ # Optim can be quite problematic, only use on small regions then find most reasonable fit over all regions
    #interceptG2ps[k]= (optim(c(0),f = function(x) (f2(x)-0.95*(f1(x)))^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    interceptG2ps[k]= (optim(c(0),f = function(x) (pnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]),lower.tail  =F)-0.025)^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    interceptG2p1s[k]= (optim(c(0),f = function(x) (pnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]),lower.tail  =F)-0.01)^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    
    if((interceptG2ps[k] > (((k+1)*0.1)-0.049))|(interceptG2ps[k] < (((k-1)*0.1)+0.049))){
      interceptG2ps[k] = 0
    }
    if((interceptG2p1s[k] > (((k+1)*0.1)-0.049))|(interceptG2p1s[k] < (((k-1)*0.1)+0.049))){
      interceptG2p1s[k] = 0
    }
  }
  alpha = 0.05/nrow(sData) # Bonferroni adjustment
  interceptG2bs =c()
  for(k in 1:10){
    interceptG2bs[k]= (optim(c(0),f = function(x) (pnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]),lower.tail  =F)-alpha)^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    if((interceptG2bs[k] > (((k+1)*0.1)-0.049))|(interceptG2bs[k] < (((k-1)*0.1)+0.049))){
      interceptG2bs[k] = 0
    }
  }
  Xp = max(interceptG2ps)
  Xp1 = max(interceptG2p1s)
  Xb = max(interceptG2bs)
  
  # Find FDR q-value cutoff
  #f3 = function(x) 1 - (f2(x)/f1(x)) # = P(Xi %in% C1 | LGR = x)
  f3 = function(x) pnorm(x, mean = modG2$parameters$mean[minG], sd = sqrt(modG2$parameters$variance$sigmasq[minG]),lower.tail  =F)
  
  qData = sData %>% mutate(qvals = f3(Best.LGR)) %>% arrange(qvals)
  
  qData = qData %>% mutate(qcomparison = (1:nrow(qData))*0.05/nrow(qData), success = qvals<=qcomparison)
  
  XpFDR = min(filter(qData, success, Best.LGR>0)$Best.LGR)
  
  fa = function(x) pros[minG]*dnorm(x, mean = means[minG], sd = sds[minG])
  fb = function(x) pros[maxG]*dnorm(x, mean = means[maxG], sd = sds[maxG])
  fc = function(x) fa(x)+fb(x)
  predictedcutoff = means[minG] + 2*sds[minG]
  
  pval = Vectorize(function(x)  1-((1-pnorm(predictedcutoff, mean = means[minG], sd = sds[minG]))*fa(x)/fc(x) + (1-pnorm(predictedcutoff, mean = x, sd = sds[maxG]))*fb(x)/fc(x)))
  interceptpV20s =c()
  interceptpV40s =c()
  for(k in 1:20){
    interceptpV20s[k]= (optim(c(0),f = function(x) (pval(x)-0.2)^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    interceptpV40s[k]= (optim(c(0),f = function(x) (pval(x)-0.4)^2, method = "Brent", lower = ((k-1)*0.1), upper = ((k+1)*0.1)))$par
    if((interceptpV20s[k] > (((k+1)*0.1)-0.049))|(interceptpV20s[k] < (((k-1)*0.1)+0.049))){
      interceptpV20s[k] = 0
    }
    if((interceptpV40s[k] > (((k+1)*0.1)-0.049))|(interceptpV40s[k] < (((k-1)*0.1)+0.049))){
      interceptpV40s[k] = 0
    }
  }
  XpV20 = max(interceptpV20s)
  XpV40 = max(interceptpV40s)
  y = arrange(sData, desc(Best.LGR))$Best.LGR
  n = 1
  x = seq(min(sData$Best.LGR),max(sData$Best.LGR),length = 1000)
  
  # Plots data
  dist = f1(x)
  if(fitPlot){
    data = hist(sData$Best.LGR, breaks = 150, plot = FALSE)
    distVals = tibble(X = x, Y = dist, Distribution = "Mixture model")
    dataVals = tibble(X = data$mids, Y = data$density, Distribution = "Screen data")
    Vals = rbind(distVals, dataVals)
   
  png(filename = paste0(Screen,"DPfit.png"), width = 600)
  print(ggplot(Vals, aes(x= X, y = Y, colour = Distribution, size  = Distribution))+ geom_line(linetype = 1) +  theme_bw() +scale_color_brewer(palette="Set1")+labs(x="LGR",y="density")+theme(text = element_text(size = 24,family = "Helvetica")) + geom_vline(aes(xintercept = Xq50,linetype = "q(x) = 50%"))+geom_vline(aes(xintercept = Zcutoff,linetype = "Z-score")) + scale_linetype_manual(name = "Cutoff", values = c("q(x) = 50%" = 2, "Z-score" = 3))+scale_size_manual(values=c(2,1))
  )
  dev.off()}

  if(compDistPlot){
    comp = tibble(X =x, Y = modG2$parameters$pro[1]*dnorm(x, mean = modG2$parameters$mean[1], sd = sqrt(modG2$parameters$variance$sigmasq[1])), Distribution = "Distribution 1")
    
    int= tibble(X=x, Y = modG2$parameters$pro[2]*dnorm(x, mean = modG2$parameters$mean[2], sd = sqrt(modG2$parameters$variance$sigmasq[2])), Distribution = paste("Distribution", 2))
    comp = rbind(comp, int)
    
    comp = rbind(comp, tibble(X=x, Y=dist, Distribution = "Mixture model"))
    
    png(filename = paste0(Screen,"DPcompdist.png"), width = 600)
    print(ggplot(comp, aes(x = X, y = Y, colour = Distribution, size  = Distribution)) + geom_line() +  theme_bw() +scale_color_brewer(palette="Set1")+labs(x="LGR",y="density")+theme(text = element_text(size = 24,family = "Helvetica"),legend.background = element_rect(fill="transparent")) + geom_vline(aes(xintercept = Xq50,linetype = "q(x) = 50%"))+geom_vline(aes(xintercept = Zcutoff,linetype = "Z-score")) + scale_linetype_manual(name = "Cutoff", values = c("q(x) = 50%" = 2, `Z-score` = 3))+scale_size_manual(values=c(1,1,1,2))
    )
    dev.off()
    
  }

  return(list(Zcutoff=Zcutoff, means=means, pros=pros, sds=sds, Xp = Xp,Xp1 = Xp1, Xb = Xb, XpFDR = XpFDR,Xq5 = Xq5,Xq50 = Xq50, XpV20 = XpV20, XpV40=XpV40,failure = failure, XpV20 = XpV20, XpV40=XpV40))
  
}

predValidation = function(ValData, Params, Screen, combineplates = F){
  # ValData should have cols ORF, LGR, Validation
  
  if("Validated?"%in%colnames(ValData)){
    ValData = rename(ValData, Validation = `Validated?`)
  }
  ValData = arrange(ValData, desc(sLGR))
  
  N = max(ValData$Plate)
  Data2 = data.frame(X = c(),Y = c(), line = c())
  BinSize = 10
  for(j in 1:N){
    plate = filter(ValData, Plate == j)
    FPR = c()
    LGR = c()
    BinSize = 10
    for(k in (BinSize+1):(nrow(plate)-BinSize)){
      int = plate[(k-BinSize):(k+BinSize),]
      FPR[k-BinSize] = nrow(filter(int,  Validation == 0))/nrow(int)
      LGR[k-BinSize] = mean(int$sLGR)
    }
    Data1 = tibble(X = LGR,Y = FPR, line = paste0("Plate", j))
    Data2 = rbind(Data2, Data1)
    
  }
  
  fa = function(x) Params$pros[1]*dnorm(x, mean = Params$means[1], sd = Params$sds[1])
  fb = function(x) Params$pros[2]*dnorm(x, mean = Params$means[2], sd = Params$sds[2])
  fc = function(x) fa(x)+fb(x)
  predictedcutoff = Params$means[1] + 2*Params$sds[1]
  
  pval = function(x)  1- ((1-pnorm(predictedcutoff, mean = Params$means[1], sd = Params$sds[1]))*fa(x)/fc(x) + (1-pnorm(predictedcutoff, mean = x, sd = Params$sds[2]))*fb(x)/fc(x))
  
  if(combineplates == TRUE){
    Data2$line = "data"
    Data2 = arrange(Data2, desc(X))
  }

  Data3 = rbind(Data2, tibble(X = seq(from = 0, to = 2, length.out = 100), Y= pval(seq(from = 0, to =2, length.out = 100)) , line = "Predicted"))
  
  
  
  png(paste0(Screen,"FPRcutoffs.png"), width = 800)
  print(ggplot(Data3, aes(x= X,y = Y, colour = line))+geom_line() + theme_bw()  +labs(x = "LGR", y= "False Positive Rate")+theme(text = element_text(size = 24,family = "Helvetica")) + ylim(c(0,1)) +xlim(c(0,2)) + geom_vline(aes(xintercept = params$XpV20, linetype = "20%"))+ geom_vline(aes(xintercept = params$XpV40, linetype = "40%")) + scale_linetype_manual(name = "FPR", values = c(2,3))
  )
  dev.off()
}

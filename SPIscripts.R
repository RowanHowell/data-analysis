AnalyseSPIData = function(GOIfile,GBPfile,plasmid,write = FALSE,writename = "SPIdata.csv", smooth = T){
  #GOIfile and GBP file should be in form "Nud1GOIcontrol.txt"
  #plasmid should be in form "pHT584"
  # writename should be eg "Nud1data.csv" - will export as .csv file
  # Note this program assumes GOIfile and GBPfile are exactly same length and order
  
  oldw <- getOption("warn") # turn off warnings for nice output!
  options(warn = -1)
#-------------------------------------------------------------------------------------------
# Section 1: load in files, merge into single file and tidy up column names and order
  
  C=c("Query","Condition","Plate #","Row","Column","P-Value","Z-Score","Normalized Growth Ratio (Comparer::Exp)","Growth Ratio (Comparer / Exp)","Log Growth Ratio","Normalized Colony Size 1","Normalized Colony Size 2","Normalized Colony Size 3","Normalized Colony Size 4","Colony Circularity 1","Colony Circularity 2","Colony Circularity 3","Colony Circularity 4","ID Column")
  
  # Fnd skip number
  GOI1 = read_delim(GOIfile,"\t", col_names = F)
  GBP1 = read_delim(GBPfile,"\t", col_names = F)
  
  GOIQs = grep("Query",GOI1$X1)
  GBPQs = grep("Query",GBP1$X1)
  
  skipGOI = GOIQs[2] # to troubleshoot try -1
  skipGBP = GBPQs[2]
  
  GOItable1 = read_delim(GOIfile,"\t",skip = skipGOI,col_types = cols(), col_names = C) #Load files into R
  GBPtable1 = read_delim(GBPfile,"\t",skip = skipGBP,col_types = cols(), col_names = C)
  
  GOItable2 = filter(GOItable1, Query == plasmid) #Select only experimental comparison
  GBPtable2 = filter(GBPtable1, Query == plasmid)
  
  # Some names cause problems, especially when they include " "
  GOItable2 = rename(GOItable2, GOI.Z = 'Z-Score', GOI.LGR = 'Log Growth Ratio', Plate = 'Plate #', ORF = 'ID Column', NMR = 'Normalized Growth Ratio (Comparer::Exp)')      
  GBPtable2 = rename(GBPtable2, GBP.Z = 'Z-Score', GBP.LGR = 'Log Growth Ratio', Plate = 'Plate #', ORF = 'ID Column', NMR = 'Normalized Growth Ratio (Comparer::Exp)')
  
  MergeFile = mutate(GOItable2, GBP.Z = GBPtable2$GBP.Z, GBP.LGR = GBPtable2$GBP.LGR) # add GBP cols
  MergeFile = select(MergeFile, Query, Plate, Row, Column, ORF, GOI.Z, GOI.LGR, GBP.Z, GBP.LGR) # remove unnecessary cols
  
  MergeFile = mutate(MergeFile, Plate = gsub("\\[|\\]", "", Plate)) # Correct plate col to numbers
#------------------------------------------------------------------------------------------- 
# Section 2: Sort out colony size column and create dead and excluded Booleans
  
  GOI.NMR = strsplit(GOItable2$NMR,"::") #Split column at "::"
  GBP.NMR = strsplit(GBPtable2$NMR,"::")
  
  GOI.NMRcomp = unlist(lapply(GOI.NMR,'[',1)) # First entries of split string list
  GOI.NMRexp = unlist(lapply(GOI.NMR,'[',2)) # Second entries of split string list
  GOI.NMRdead = numeric(length(GOI.NMRexp))
  GOI.NMRexcl = GOI.NMRdead
  
  GBP.NMRcomp = unlist(lapply(GBP.NMR,'[',1)) # First entries of split string list
  GBP.NMRexp = unlist(lapply(GBP.NMR,'[',2)) # Second entries of split string list
  GBP.NMRdead = numeric(length(GOI.NMRexp))
  GBP.NMRexcl = GOI.NMRdead
  
  for(i in 1:length(GOI.NMRexp)){
    if(substr(GOI.NMRcomp[i],1,5) == "dead-"){
      GOI.NMRcomp[i] = substr(GOI.NMRcomp[i],6,9) #remove "dead-" part
      GOI.NMRdead[i] = 1 #flag as dead colony
    }
    else if(substr(GOI.NMRcomp[i],1,9) == "excluded-"){
      GOI.NMRcomp[i] = substr(GOI.NMRcomp[i],10,13) #remove "excluded-" part
      GOI.NMRexcl[i] = 1 #flag as excluded colony
    }
    
    if(substr(GOI.NMRexp[i],5,9) == "-dead"){
      GOI.NMRexp[i] = substr(GOI.NMRexp[i],1,4) #remove "-dead" part
    }
    else if(substr(GOI.NMRexp[i],5,13) == "-excluded"){
      GOI.NMRexp[i] = substr(GOI.NMRexp[i],1,4) #remove "excluded-" part
    }
    
    
    if(substr(GBP.NMRcomp[i],1,5) == "dead-"){
      GBP.NMRcomp[i] = substr(GBP.NMRcomp[i],6,9) #remove "dead-" part
      GBP.NMRdead[i] = 1 #flag as dead colony
    }
    else if(substr(GBP.NMRcomp[i],1,9) == "excluded-"){
      GBP.NMRcomp[i] = substr(GBP.NMRcomp[i],10,13) #remove "excluded-" part
      GBP.NMRexcl[i] = 1 #flag as excluded colony
    }
    
    if(substr(GBP.NMRexp[i],5,9) == "-dead"){
      GBP.NMRexp[i] = substr(GBP.NMRexp[i],1,4) #remove "-dead" part
    }
    else if(substr(GBP.NMRexp[i],5,13) == "-excluded"){
      GBP.NMRexp[i] = substr(GBP.NMRexp[i],1,4) #remove "excluded-" part
    }
  }
  GOI.NMRcomp = as.numeric(GOI.NMRcomp) # Force R to treat this as double col
  GOI.NMRexp = as.numeric(GOI.NMRexp)
  GBP.NMRcomp = as.numeric(GBP.NMRcomp)
  GBP.NMRexp = as.numeric(GBP.NMRexp)
  
  MergeFile = mutate(MergeFile,GOI.NMRcomp,GOI.NMRexp,GOI.NMRdead,GOI.NMRexcl,GBP.NMRcomp,GBP.NMRexp,GBP.NMRdead,GBP.NMRexcl) 
  
#-------------------------------------------------------------------------------------------
# Section 3: Clean data, cut blank/failed columns
  
  MergeFile = filter(MergeFile, (GOI.NMRdead+GOI.NMRexcl == 0)|(GBP.NMRdead+GBP.NMRexcl == 0) ) # Ensure either GOI OR GBP comparison has good data
  
  MergeFile = filter(MergeFile, !(ORF %in% c("BLANK","Unknown1","Unknown2"))) # Remove Blank spaces and Unknowns
  MergeFile = arrange(MergeFile,ORF)
#-------------------------------------------------------------------------------------------  
# Section 4: Smooth data
  if (smooth) {
    MergeFileS = SmoothSPIdata(MergeFile)
    MergeFileS = arrange(MergeFileS,ORF)
    MergeFile = MergeFile %>% rename(usGOI.Z = GOI.Z, usGBP.Z = GBP.Z) %>% add_column(GBP.Z = MergeFileS$`GBP.Z corrected scores`,GOI.Z = MergeFileS$`GOI.Z corrected scores`)
  }
  
#-------------------------------------------------------------------------------------------  
# Section 5: Make best Z-score/LGR columns
  
  # Make average Z-score column, if for one comparison, data is n/a use the other comparison
  B1 = as.numeric(!(MergeFile$GOI.NMRdead+MergeFile$GOI.NMRexcl)) # B1[j] = 1 when jth entry is good data, 0 otherwise
  B2 = as.numeric(!(MergeFile$GBP.NMRdead+MergeFile$GBP.NMRexcl)) 
  B3 = B1*B2 #B3[j] = 1 when both B1 and B2 are = 1, 0 otherwise
  u = (as.numeric(MergeFile$GOI.Z)+as.numeric(MergeFile$GBP.Z))/2 # mean of columns
  MergeFile = mutate(MergeFile,Best.Z = ifelse(B3,u, ifelse(B1,as.numeric(MergeFile$GOI.Z),as.numeric(MergeFile$GBP.Z))))
  
  # Make average LGR column
  v = (as.numeric(MergeFile$GOI.LGR)+as.numeric(MergeFile$GBP.LGR))/2 # mean of columns
  MergeFile = mutate(MergeFile,Best.LGR = ifelse(B3,v, ifelse(B1,as.numeric(MergeFile$GOI.LGR),as.numeric(MergeFile$GBP.LGR))))

  MergeFile = MergeFile %>% arrange(Plate, Row, Column)
#-------------------------------------------------------------------------------------------
# Section 6: write data to file
  
  if(write == TRUE){ # Write data to file if desired
    write_csv(MergeFile,writename)
  }
  options(warn = oldw) # Switch warnings back on
  
  return(MergeFile)
}

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

SPIplots = function(MergeFile){
  MergeFile = MergeFile %>% filter(!(is.na(GBP.Z)))
  A = arrange(MergeFile, desc(Best.Z)) #arrange largest to smallest
  lim = ceiling(max(A$Best.Z)) # find best limits for graph
  plot(A$Best.Z, xlab = "GFP strains", ylab = " Z-score", pch = 19,col = rainbow(8000),cex = .4,panel.first=grid(20,10),ylim = c(-1*lim,lim))
  abline(h = c(2,-2)) # add lines at ±2
  
  
  B = filter(MergeFile, !(is.na(GOI.Z))&!(is.na(GBP.Z)))
  Xlim = ceiling(as.numeric(max(B$GOI.Z))) # Best limits for graph
  Ylim = ceiling(as.numeric(max(B$GBP.Z)))
  Color = 9*(B$Best.Z>=2)+1 # Vector of colors, red = SPI, black  = non-SPI
  
  plot(B$GOI.Z,B$GBP.Z,pch = 19,panel.first=grid(12,12),asp =1,cex = .4,ylim = c(-1*Ylim,Ylim),xlim = c(-1*Xlim,Xlim), xlab = "GOI control Z-score", ylab = "GBP control Z-score", col = Color, main = "Control Correlation" )
  lines(c(100,-96),c(-96,100)) # Line showing (x+y)/2 = 2, line showing SPI limit
  abline(h = 2,v=2)
  
  legend("bottomleft",c("SPI","non-SPI"), col = c(10,1), pch =19)
}

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

SmoothSPIdata = function(data,output){
  
  data = data %>% select(plate = Plate, row = Row, column = Column, ORF, GBP.Z, GOI.Z)
  write_tsv(data, "results.txt")
  system(paste("perl", "2ndRCsmoother-16.plx"))
  
  smoothed = read_tsv("smoothed.tab")
  return(smoothed)
}

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

CheckSmoothing = function(smoothed){
  
  data = arrange(smoothed, Plate, Row, as.numeric(Column))
  
  par(mfrow=c(4, 1))
  plot(smoothed$usGOI.Z, pch = 19, col = "Red",cex = .2)
  plot(smoothed$GOI.Z, pch = 19, col = "Black",cex = .2)
  plot(smoothed$usGBP.Z, pch = 19, col = "Red",cex = .2)
  plot(smoothed$GBP.Z, pch = 19, col = "Black",cex = .2)
  
}

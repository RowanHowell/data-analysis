AnalyseSPIData = function(GOIfile,GBPfile,plasmid,write = FALSE,writename = "SPIdata.csv", smooth = T, smoothLGR = F, aggregate =FALSE){
  # GOIfile and GBP file should be in form "Nud1GOIcontrol.txt"
  # plasmid should be in form "pHT584"
  # writename should be eg "Nud1data.csv" - will export as .csv file
  # smooth = TRUE will smooth Z-scores, you need Peter's 2ndRCsmoother-16.plx file in the same folder!
  # smoothLGR = TRUE will smooth LGRs too
  # aggregate = TRUE will output a reduced file where duplicate colonies are averaged and the only cols are ORF, Best.Z and Best.LGR
  # Latest version 6/9/18, fixed smoothing issues 
  # Latest version 14/9/18, add LGR smoothing 
  # Latest version 24/9/18, fixed gap issues, added aggregate option, fixed some bugs 
  
  
  oldw <- getOption("warn") # turn off warnings for nice output!
  options(warn = -1)
  gap = 1 #Try gap =1 first
#-------------------------------------------------------------------------------------------
# Section 1: load in files, merge into single file and tidy up column names and order
  
  # Fnd skip number
  GOI1 = read_delim(GOIfile,"\t", col_names = F)
  GBP1 = read_delim(GBPfile,"\t", col_names = F)
  
  GOIQs = grep("Query",GOI1$X1)
  GBPQs = grep("Query",GBP1$X1)
  
  skipGOI = max(GOIQs)-gap
  skipGBP = max(GBPQs)-gap
  
  GOItable1 = read_delim(GOIfile,"\t",skip = skipGOI,col_types = cols()) #Load files into R
  GBPtable1 = read_delim(GBPfile,"\t",skip = skipGBP,col_types = cols())
  
  if(colnames(GOItable1)[1]=="Query"){ # Check if properly read
    GOItable2 = filter(GOItable1, Query == plasmid) #Select only experimental comparison
  } else{
    gap2 = abs(gap-1) # switch 0 to 1 or 1 to 0
    skipGOI = max(GOIQs)-gap2
    GOItable1 = read_delim(GOIfile,"\t",skip = skipGOI,col_types = cols()) 
    GOItable2 = filter(GOItable1, Query == plasmid)
  }
  
  if(colnames(GBPtable1)[1]=="Query"){
    GBPtable2 = filter(GBPtable1, Query == plasmid) #Select only experimental comparison
  } else{
    gap2 = abs(gap-1) # switch 0 to 1 or 1 to 0
    skipGBP = max(GBPQs)-gap2
    GBPtable1 = read_delim(GBPfile,"\t",skip = skipGBP,col_types = cols()) #Load files into R
    GBPtable2 = filter(GBPtable1, Query == plasmid)
  }
  
  
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
  GBP.NMRdead = numeric(length(GBP.NMRexp))
  GBP.NMRexcl = GBP.NMRdead
  
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
    MergeFileS = select(MergeFileS,Plate, Row, Column, ORF, GBP.Z = `GBP.Z corrected scores`, GOI.Z = `GOI.Z corrected scores`) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    MergeFile = MergeFile %>% rename(usGOI.Z = GOI.Z, usGBP.Z = GBP.Z) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    
    MergeFile = inner_join(MergeFile, MergeFileS, by = c("Plate", "Row", "Column", "ORF"))
  }
  
# Section 4.1 : Smooth LGRdata
  if (smoothLGR) {
    MergeFileS2 = SmoothSPIdata(MergeFile, Z = FALSE)
    MergeFileS2 = select(MergeFileS2,Plate, Row, Column, ORF, GBP.LGR = `GBP.LGR corrected scores`, GOI.LGR = `GOI.LGR corrected scores`) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    MergeFile = MergeFile %>% rename(usGOI.LGR = GOI.LGR, usGBP.LGR = GBP.LGR) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    
    MergeFile = inner_join(MergeFile, MergeFileS2, by = c("Plate", "Row", "Column", "ORF"))
  }
  
#-------------------------------------------------------------------------------------------  
# Section 5: Make best Z-score/LGR columns
  
  # Make average Z-score column, if for one comparison, data is n/a use the other comparison
  B1 = as.numeric(!(MergeFile$GOI.NMRdead+MergeFile$GOI.NMRexcl)) # B1[j] = 1 when jth entry is good data, 0 otherwise
  B2 = as.numeric(!(MergeFile$GBP.NMRdead+MergeFile$GBP.NMRexcl)) 
  B3 = B1*B2 #B3[j] = 1 when both B1 and B2 are = 1, 0 otherwise
  u = (as.numeric(MergeFile$GOI.Z)+as.numeric(MergeFile$GBP.Z))/2 # mean of columns
  MergeFile = mutate(MergeFile,Best.Z = ifelse(B3,u, ifelse(B1,as.numeric(MergeFile$GOI.Z),as.numeric(MergeFile$GBP.Z))))
  
  # Make average Z-score column for us data, if for one comparison, data is n/a use the other comparison
  if(smooth){
    u = (as.numeric(MergeFile$usGOI.Z)+as.numeric(MergeFile$usGBP.Z))/2 # mean of columns
    MergeFile = mutate(MergeFile,usBest.Z = ifelse(B3,u, ifelse(B1,as.numeric(MergeFile$usGOI.Z),as.numeric(MergeFile$usGBP.Z))))
  }
  # Make average LGR column
  v = (as.numeric(MergeFile$GOI.LGR)+as.numeric(MergeFile$GBP.LGR))/2 # mean of columns
  MergeFile = mutate(MergeFile,Best.LGR = ifelse(B3,v, ifelse(B1,as.numeric(MergeFile$GOI.LGR),as.numeric(MergeFile$GBP.LGR))))
  
  # Make average LGR column for us data
  if(smoothLGR){
    v = (as.numeric(MergeFile$usGOI.LGR)+as.numeric(MergeFile$usGBP.LGR))/2 # mean of columns
    MergeFile = mutate(MergeFile,usBest.LGR = ifelse(B3,v, ifelse(B1,as.numeric(MergeFile$usGOI.LGR),as.numeric(MergeFile$usGBP.LGR))))
  }
  
  MergeFile = MergeFile %>% arrange(Plate, Row, Column)
#-------------------------------------------------------------------------------------------
# Section 6: Aggregate data
# If aggregate = TRUE a reduced version of the analysis will be returned with unique ORFs
  if(aggregate){
    LGRData = aggregate(MergeFile$Best.LGR,by=list(ORF=MergeFile$ORF),FUN=mean) %>% rename(Best.LGR = x)
    ZData = aggregate(MergeFile$Best.Z,by=list(ORF=MergeFile$ORF),FUN=mean)%>% rename(Best.Z = x)
  
    MergeFile = inner_join(ZData, LGRData, by = "ORF")
  }
#-------------------------------------------------------------------------------------------
# Section 7: write data to file
  
  if(write == TRUE){ # Write data to file if desired
    write_csv(MergeFile,writename)
  }
  options(warn = oldw) # Switch warnings back on
  
  return(MergeFile)
}

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

SmoothSPIdata = function(data,output, Z = TRUE, ctrl = FALSE){
  if(!ctrl){
    if(Z){ # Smoothes Z by default
      data = data %>% dplyr::select(plate = Plate, row = Row, column = Column, ORF, GBP.Z, GOI.Z)
    }
    else{ # Smoothes LGRs
      data = data %>% dplyr::select(plate = Plate, row = Row, column = Column, ORF, GBP.LGR, GOI.LGR)
    }
  }else{
    if(Z){ # Smoothes Z by default
      data = data %>% dplyr::select(plate = Plate, row = Row, column = Column, ORF, GOI.Z)
    }
    else{ # Smoothes LGRs
      data = data %>% dplyr::select(plate = Plate, row = Row, column = Column, ORF, GOI.LGR)
    }
  }
  write_tsv(data, "results.txt")
  system(paste("perl", "2ndRCsmoother-16.plx"))
  
  smoothed = read_tsv("smoothed.tab")
  return(smoothed)
}
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

AnalyseValidationData = function(GOIfile,GBPfile,plasmid,write = FALSE,writename = "SPIdata.csv",K=1){
    # This is a version of AnalyseSPIdata adapted for validation screens
    # Key changes are that this accounts for some variations in plasmid naming
    # It can also be used to synthesise validation experiments from multiple plates
    # Some features from AnalyseSPIdata like aggregation and smoothing have been removed as they are not appropriate to validation
    # GOIfile and GBP file should be in form "Nud1GOIcontrol.txt"
    # plasmid should be in form "pHT584"
    # writename should be eg "Nud1data.csv" - will export as .csv file
    # Latest version 6/9/18, fixed smoothing issues 
    # Latest version 14/9/18, add LGR smoothing 
    # Latest version 24/9/18, fixed gap issues, added aggregate option, fixed some bugs 
    
    
    oldw <- getOption("warn") # turn off warnings for nice output!
    options(warn = -1)
    gap = 0 #Try gap =1 first
    #-------------------------------------------------------------------------------------------
    # Section 1: load in files, merge into single file and tidy up column names and order
    
    # Fnd skip number
    GOI1 = read_delim(GOIfile,"\t", col_names = F)
    GBP1 = read_delim(GBPfile,"\t", col_names = F)
    
    GOIQs = grep("Query",GOI1$X1)
    GBPQs = grep("Query",GBP1$X1)
    
    skipGOI = max(GOIQs)-gap
    skipGBP = max(GBPQs)-gap
    
    GOItable1 = read_delim(GOIfile,"\t",skip = skipGOI,col_types = cols()) #Load files into R
    GBPtable1 = read_delim(GBPfile,"\t",skip = skipGBP,col_types = cols())
    
    #Note next bit deviates from main algorithm
    if(colnames(GOItable1)[1]!="Query"){ # Check if properly read
      gap2 = abs(gap-1) # switch 0 to 1 or 1 to 0
      skipGOI = max(GOIQs)-gap2
      GOItable1 = read_delim(GOIfile,"\t",skip = skipGOI,col_types = cols())
    }
    vals = unique(GOItable1$Query) # This is necessary to account for variations in naming
    if(!(plasmid %in% vals)){
      plas = substr(plasmid, 4, 6)
      if(!(plas %in% vals)){
        plas = paste0(plas,plateID[K])
      }
      plasmid = plas
    }
    
    GOItable2 = filter(GOItable1, Query == plasmid)
    
    
    if(colnames(GBPtable1)[1]!="Query"){ # Check if properly read
      gap2 = abs(gap-1) # switch 0 to 1 or 1 to 0
      skipGBP = max(GBPQs)-gap2
      GBPtable1 = read_delim(GBPfile,"\t",skip = skipGBP,col_types = cols())
    }
    GBPtable2 = filter(GBPtable1, Query == plasmid)
    
    
    
    # Some names cause problems, especially when they include " "
    GOItable2 = rename(GOItable2, GOI.LGR = 'Log Growth Ratio', Plate = 'Plate #', ORF = 'ID Column', NMR = 'Normalized Growth Ratio (Comparer::Exp)')      
    GBPtable2 = rename(GBPtable2, GBP.LGR = 'Log Growth Ratio', Plate = 'Plate #', ORF = 'ID Column', NMR = 'Normalized Growth Ratio (Comparer::Exp)')
    
    MergeFile = mutate(GOItable2, GBP.LGR = GBPtable2$GBP.LGR) # add GBP cols
    MergeFile = select(MergeFile, Query, Plate, Row, Column, ORF,GOI.LGR,GBP.LGR) # remove unnecessary cols
    
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
    GBP.NMRdead = numeric(length(GBP.NMRexp))
    GBP.NMRexcl = GBP.NMRdead
    
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
    # Section 4: Make best Z-score/LGR columns
    
    # Make average Z-score column, if for one comparison, data is n/a use the other comparison
    B1 = as.numeric(!(MergeFile$GOI.NMRdead+MergeFile$GOI.NMRexcl)) # B1[j] = 1 when jth entry is good data, 0 otherwise
    B2 = as.numeric(!(MergeFile$GBP.NMRdead+MergeFile$GBP.NMRexcl)) 
    B3 = B1*B2 #B3[j] = 1 when both B1 and B2 are = 1, 0 otherwise
    
    # Make average LGR column
    v = (as.numeric(MergeFile$GOI.LGR)+as.numeric(MergeFile$GBP.LGR))/2 # mean of columns
    MergeFile = mutate(MergeFile,Best.LGR = ifelse(B3,v, ifelse(B1,as.numeric(MergeFile$GOI.LGR),as.numeric(MergeFile$GBP.LGR))))
    
    MergeFile = MergeFile %>% arrange(Plate, Row, Column)
    #-------------------------------------------------------------------------------------------
    # Section 5: write data to file
    
    if(write == TRUE){ # Write data to file if desired
      write_csv(MergeFile,writename)
    }
    options(warn = oldw) # Switch warnings back on
    
    return(MergeFile)
}

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

AnalyseCTRLData = function(GOIfile,write = FALSE,writename = "SPICTRLdata.csv", smooth = T, smoothLGR = F, aggregate =FALSE, plasmid = "pHT4"){
  # This version of AnalyseSPIdata can be used to extract "null comparisons" eg LGR calculated by comparing controls
  # This is VERY UNLIKELY to be the program you want unless you know what you are doing.
  # GOIfile and GBP file should be in form "Nud1GOIcontrol.txt"
  # plasmid should be in form "pHT584"
  # writename should be eg "Nud1data.csv" - will export as .csv file
  # Latest version 6/9/18, fixed smoothing issues 
  # Latest version 14/9/18, add LGR smoothing 
  # Latest version 24/9/18, fixed gap issues, added aggregate option, fixed some bugs 
  
  select = dplyr::select
  
  oldw <- getOption("warn") # turn off warnings for nice output!
  options(warn = -1)
  gap = 1 #Try gap =1 first
  #-------------------------------------------------------------------------------------------
  # Section 1: load in files, merge into single file and tidy up column names and order
  
  # Fnd skip number
  GOI1 = read_delim(GOIfile,"\t", col_names = F)

  GOIQs = grep("Query",GOI1$X1)

  skipGOI = max(GOIQs)-gap

  GOItable1 = read_delim(GOIfile,"\t",skip = skipGOI,col_types = cols()) #Load files into R

  if(colnames(GOItable1)[1]=="Query"){ # Check if properly read
    GOItable2 = filter(GOItable1, Query == plasmid) #Select only experimental comparison
  } else{
    gap2 = abs(gap-1) # switch 0 to 1 or 1 to 0
    skipGOI = max(GOIQs)-gap2
    GOItable1 = read_delim(GOIfile,"\t",skip = skipGOI,col_types = cols()) 
    GOItable2 = filter(GOItable1, Query == plasmid)
  }
  
  if("Normalized Ratio (Comparer::Exp)" %in% colnames(GOItable2)){
    GOItable2 = rename(GOItable2,`Normalized Growth Ratio (Comparer::Exp)` = `Normalized Ratio (Comparer::Exp)`)
  }
  if("Calculated Log Ratio (Comparer::Exp)" %in% colnames(GOItable2)){
    GOItable2 = rename(GOItable2,`Log Growth Ratio` = `Calculated Log Ratio (Comparer::Exp)`)
  }
  if("Col" %in% colnames(GOItable2)){
    GOItable2 = rename(GOItable2,`Column` = `Col`)
  }

  
  # Some names cause problems, especially when they include " "
  GOItable2 = rename(GOItable2, GOI.Z = 'Z-Score', GOI.LGR = 'Log Growth Ratio', Plate = 'Plate #', ORF = 'ID Column', NMR = 'Normalized Growth Ratio (Comparer::Exp)')      

  MergeFile = GOItable2 
  MergeFile = dplyr::select(MergeFile, Query, Plate, Row, Column, ORF, GOI.Z, GOI.LGR) # remove unnecessary cols
  MergeFile = mutate(MergeFile, Plate = gsub("\\[|\\]", "", Plate)) # Correct plate col to numbers
  #------------------------------------------------------------------------------------------- 
  # Section 2: Sort out colony size column and create dead and excluded Booleans
  
  GOI.NMR = strsplit(GOItable2$NMR,"::") #Split column at "::"

  GOI.NMRcomp = unlist(lapply(GOI.NMR,'[',1)) # First entries of split string list
  GOI.NMRexp = unlist(lapply(GOI.NMR,'[',2)) # Second entries of split string list
  GOI.NMRdead = numeric(length(GOI.NMRexp))
  GOI.NMRexcl = GOI.NMRdead
  
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
    
  }
  GOI.NMRcomp = as.numeric(GOI.NMRcomp) # Force R to treat this as double col
  GOI.NMRexp = as.numeric(GOI.NMRexp)

  MergeFile = mutate(MergeFile,GOI.NMRcomp,GOI.NMRexp,GOI.NMRdead,GOI.NMRexcl) 
  
  #-------------------------------------------------------------------------------------------
  # Section 3: Clean data, cut blank/failed columns
  
  MergeFile = filter(MergeFile, (GOI.NMRdead+GOI.NMRexcl == 0)) # Ensure either GOI OR GBP comparison has good data
  
  MergeFile = filter(MergeFile, !(ORF %in% c("BLANK","Unknown1","Unknown2"))) # Remove Blank spaces and Unknowns
  MergeFile = arrange(MergeFile,ORF)
  #-------------------------------------------------------------------------------------------  
  # Section 4: Smooth data
  if (smooth) {
    MergeFileS = SmoothSPIdata(MergeFile, ctrl = TRUE)
    MergeFileS = select(MergeFileS,Plate, Row, Column, ORF, GOI.Z = `GOI.Z corrected scores`) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    MergeFile = MergeFile %>% rename(usGOI.Z = GOI.Z) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    
    MergeFile = inner_join(MergeFile, MergeFileS, by = c("Plate", "Row", "Column", "ORF"))
  }
  
  # Section 4.1 : Smooth LGRdata
  if (smoothLGR) {
    MergeFileS2 = SmoothSPIdata(MergeFile, Z = FALSE, ctrl = TRUE)
    MergeFileS2 = select(MergeFileS2,Plate, Row, Column, ORF, GOI.LGR = `GOI.LGR corrected scores`) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    MergeFile = MergeFile %>% rename(usGOI.LGR = GOI.LGR) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
    
    MergeFile = inner_join(MergeFile, MergeFileS2, by = c("Plate", "Row", "Column", "ORF"))
  }
  
  MergeFile = rename(MergeFile, Best.Z = GOI.Z, Best.LGR = GOI.LGR)
  #-------------------------------------------------------------------------------------------
  # Section 6: Aggregate data
  # If aggregate = TRUE a reduced version of the analysis will be returned with unique ORFs
  if(aggregate){
    LGRData = aggregate(MergeFile$Best.LGR,by=list(ORF=MergeFile$ORF),FUN=mean) %>% rename(Best.LGR = x)
    ZData = aggregate(MergeFile$Best.Z,by=list(ORF=MergeFile$ORF),FUN=mean)%>% rename(Best.Z = x)
    
    MergeFile = inner_join(ZData, LGRData, by = "ORF")
  }
  #-------------------------------------------------------------------------------------------
  # Section 7: write data to file
  
  if(write == TRUE){ # Write data to file if desired
    write_csv(MergeFile,writename)
  }
  options(warn = oldw) # Switch warnings back on
  
  return(MergeFile)
}

  
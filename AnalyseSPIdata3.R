AnalyseSPIData = function(GOIfile,GBPfile,plasmid,write = FALSE,writename = "SPIdata.csv"){
  #GOIfile and GBP file should be in form "Nud1GOIcontrol.txt"
  #plasmid should be in form "pHT584"
  # writename should be eg "Nud1data.csv" - will export as .csv file
  # Note this program assumes GOIfile and GBPfile are exactly same length and order
  
  oldw <- getOption("warn") # turn off warnings for nice output!
  options(warn = -1)
#-------------------------------------------------------------------------------------------
# Section 1: load in files, merge into single file and tidy up column names and order
  
  GOItable1 = read_delim(GOIfile,"\t",skip = 16,col_types = cols()) #Load files into R
  GBPtable1 = read_delim(GBPfile,"\t",skip = 16,col_types = cols())
  
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
# Section 3: Clean data, cut blank/failed columns and make best mean column
  
  MergeFile = filter(MergeFile, (GOI.NMRdead+GOI.NMRexcl == 0)|(GBP.NMRdead+GBP.NMRexcl == 0) ) # Ensure either GOI OR GBP comparison has good data
  
  MergeFile = filter(MergeFile, ORF!="BLANK") # Remove Blank spaces
  
  # Make average Z-score column, if for one comparison, data is n/a use the other comparison
  B1 = as.numeric(!(MergeFile$GOI.NMRdead+MergeFile$GOI.NMRexcl)) # B1[j] = 1 when jth entry is good data, 0 otherwise
  B2 = as.numeric(!(MergeFile$GBP.NMRdead+MergeFile$GBP.NMRexcl)) 
  B3 = B1*B2 #B3[j] = 1 when both B1 and B2 are = 1, 0 otherwise
  u = (as.numeric(MergeFile$GOI.Z)+as.numeric(MergeFile$GBP.Z))/2 # mean of columns
  MergeFile = mutate(MergeFile,Best.Z = ifelse(B3,u, ifelse(B1,as.numeric(MergeFile$GOI.Z),as.numeric(MergeFile$GBP.Z))))
  
  # Make average LGR column
  v = (as.numeric(MergeFile$GOI.LGR)+as.numeric(MergeFile$GBP.LGR))/2 # mean of columns
  MergeFile = mutate(MergeFile,Best.LGR = ifelse(B3,v, ifelse(B1,as.numeric(MergeFile$GOI.LGR),as.numeric(MergeFile$GBP.LGR))))

#-------------------------------------------------------------------------------------------
# Section 4: write data to file
  
  if(write == TRUE){ # Write data to file if desired
    write_csv(MergeFile,writename)
  }
  options(warn = oldw) # Switch warnings back on
  
  return(MergeFile)
  
}


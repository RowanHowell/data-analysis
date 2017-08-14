# install.packages("tidyverse") make sure to install before first use
library(tidyverse)

source("AnalyseSPIData3.R") # Load functions into R
source("SPIplots.R")

GOIfile = "Abc1GOIcontrol.txt" # Ensure name is in quotes and has file extension
GBPfile = "Abc1GBPcontrol.txt"
plasmid = "pHT584"
writename = "Nud1data.csv" # This must end in .csv!

H = AnalyseSPIData(GOIfile,GBPfile,plasmid,write = TRUE, writename)
# Your data is now stored in the data frame H and has been written to file

SPIplots(H)
# THis will make a plot showing correlation of controls and spread of data
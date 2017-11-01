# install.packages("tidyverse") make sure to install before first use
library(tidyverse)

source("SPIscripts.R") # Load functions into R


GOIfile = "Abc1GOIcontrol.txt" # Ensure name is in quotes and has file extension
GBPfile = "Abc1GBPcontrol.txt"
plasmid = "pHT999"
writename = "Abc1data.csv" # This must end in .csv!

H = AnalyseSPIData(GOIfile,GBPfile,plasmid,write = TRUE, writename, smooth = T)
# Your data is now stored in the data frame H and has been written to file

SPIplots(H)
# THis will make a plot showing correlation of controls and spread of data



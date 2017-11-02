# install.packages("tidyverse") make sure to install before first use
library(tidyverse)

source("SPIscripts.R") # Load functions into R


GOIfile = "AbcGOIcontrol.txt" # Ensure name is in quotes and has file extension
GBPfile = "AbcGBPcontrol.txt"
plasmid = "pHT999"
writename = "Abc1data.csv" # This must end in .csv!

H = AnalyseSPIData(GOIfile,GBPfile,plasmid,write = TRUE, writename, smooth = TRUE)
# Your data is now stored in the data frame H and has been written to file

SPIplots(H)
# THis will make a plot showing correlation of controls and spread of data

CheckSmoothing(H)
# This will make a plot showing the plate effects for each of the comparisons before and after smoothing
# Note you may need to increase the height of the plot window in Rstudio if you see " Error in plot.new() : figure margins too large "

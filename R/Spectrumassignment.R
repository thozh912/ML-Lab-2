library(XLConnect)
#FROM THIS FILE LOCATION, EXCEL FILES SHOULD BE FOUND IN A SUBFOLDER IN THIS FILE LOCATION CALLED DATA
wb = loadWorkbook("C:/Users/Dator/Documents/R_HW/ML-Lab-2/data/NIRSpectra.xls")
data2 = readWorksheet(wb, sheet = "NIRSpectra", header = TRUE)
head(data2)
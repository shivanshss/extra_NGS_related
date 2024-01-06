#install.packages("devtools")
library(devtools)

devtools::install_github("Farre-lab/syntenyPlotteR")
library(syntenyPlotteR)

reformat.syntenyData("/home/bloodmark/workarea/synteny/LG2/sibeliaz_out/alignment.maf", "reformat_lg2_maf", directory = "/home/bloodmark/workarea/synteny/LG2/")

#Evolution Highway

draw.eh("output", chrRange, "data_file", fileformat = "png", colour = "lightblue", inverted.colour = "lightpink", w = 5.5, h = 10, ps = 10)

draw.ideogram("file_data", "sizefile", "output", fileformat = "png", colours = colours.default, w=8.5, h=10, ps=5)

draw.linear(output, sizefile, ..., fileformat = "png", colours = colours.default, w=13, h=5, opacity = .5)













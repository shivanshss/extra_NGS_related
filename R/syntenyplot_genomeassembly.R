devtools::install_github("marta-fb/syntenyPlotteR")
install.packages("tidyverse")
install.packages("stringr")

library(syntenyPlotteR)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(stringr)


update.packages(ask = FALSE)

getwd()

setwd("/home/bloodmark/Downloads/syntenyPlotteR/data")

draw.pairwise("trial_w_3","lengths.txt","alignment_1.txt","alignment_2.txt", fileformat = "png", w=13, h=5)


draw.pairwise("trial_w_2","lengths.txt","alignment_1.txt", fileformat = "png", w=13, h=5)

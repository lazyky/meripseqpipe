setwd("E:/zky/m6Apipe/viewer/")
source("ui.R")
source("server.R")
m6apipe.results.file = "E:/zky/m6Aviewer/MATK_arranged_results_20190619.m6APipe"
load(m6apipe.results.file)
runApp(appDir = "E:/zky/m6Apipe/viewer/",port = 8848,launch.browser = T)

#!/bin/Rscript
## Rscript get_htseq_matrix.R aligner_tools designfile gtf eg. Rscript get_htseq_matrix.R tophat2 designfile_single.txt
## designfile: filename, control_or_treated, input_or_ip, group(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)

args<-commandArgs(T)
m6APipe.result <- args[1]
m6AReport.dir <- args[2]

source(paste0(m6AReport.dir,"ui.R"))
source(paste0(m6AReport.dir,"server.R"))
load(m6APipe.result)

runApp(appDir = m6AReport.dir, port = 8848, launch.browser = F)

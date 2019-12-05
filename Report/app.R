#!/bin/Rscript
## Rscript get_htseq_matrix.R aligner_tools designfile gtf eg. Rscript get_htseq_matrix.R tophat2 designfile_single.txt
## designfile: filename, control_or_treated, input_or_ip, group(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)

args <- commandArgs(T)
args <- c("macs2_QNB_edgeR_arranged_results_2019-11-11.m6APipe",".")
m6APipe.result <- args[1]
m6AReport.dir <- args[2]
load("macs2_QNB_edgeR_arranged_results_2019-11-11.m6APipe")
## switch work dir for Shinyjs
setwd(m6AReport.dir)
source("server.R")
source("ui.R")
source("functions.R")



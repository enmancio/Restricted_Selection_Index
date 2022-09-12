#!/bin/bash
rm nohup.out 
nohup Rscript -e 'library(rmarkdown); rmarkdown::render("Supplementary_1.Rmd", "pdf_document")' &

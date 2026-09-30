FROM bioconductor/bioconductor_docker:devel

MAINTAINER pl.germain@gmail.com

WORKDIR /home/build/package

COPY . /home/build/package 

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN R -e "install.packages('remotes'); pkgs <- c('cowplot','ggseqlogo','ensembldb','AnnotationFilter','AnnotationHub','digest','DT','fst','Seqinfo','GenomicFeatures','ggplot2','htmlwidgets','Matrix','plotly','rtracklayer','shiny','shinycssloaders','shinydashboard','waiter','rintrojs', 'txdbmaker'); BiocManager::install(pkgs);"

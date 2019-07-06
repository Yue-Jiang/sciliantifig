FROM rocker/tidyverse:3.5.1

MAINTAINER Yue Jiang <rivehill@gmail.com>

# System dependencies for required R packages
RUN  rm -f /var/lib/dpkg/available \
  && rm -rf  /var/cache/apt/* \
  && apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    git

RUN apt-get update -qq && apt-get install -y --no-install-recommends gfortran

RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('HMM','depmixS4','dplyr','tidyr','purrr','ggplot2','ggrepel','gridExtra','scales','BAS','parallel','corrplot','pROC','exactRankTests','lmerTest','grid','R.utils'))"

RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('karyoploteR'))"

RUN mkdir /liftrroot/
WORKDIR /liftrroot/

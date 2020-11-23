FROM rocker/shiny-verse:3.6.1
 
# system libraries of general use
 
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libjpeg-dev \
    libxml2 \

# Install R Dependencies - installs the r packages you need - if this step fails you’re likely 
# missing system libraries that a package requires
 
 
# copy shiny-server.sh to image
 
COPY shiny-server.sh /usr/bin/
  
# copy shiny server config to image
 
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
 
# copy installer

copy installer.R /srv/shiny-server/

# copy the contents of app folder to image
 
COPY R /srv/shiny-server/R/

# install R packages

RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_version('shiny', version = '1.4.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinyBS', version = '0.61', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('mice', version = '3.8.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('DT', version = '0.12', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('ggplot2', version = '3.2.1', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('data.table', version = '1.12.8', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('truncdist', version = '1.0-2', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinyalert', version = '1.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinydashboard', version = '0.7.1', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinycssloaders', version = '0.3', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('BiocManager', version = '1.30.10', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinyjs', version = '2.0.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('shinyWidgets', version = '0.5.4', repos = 'http://cran.us.r-project.org')"
# RUN R -e "devtools::install_version('jpeg', version = '0.1-8.1', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('rexposome')"
RUN R -e "BiocManager::install('omicRexposome')"
RUN R -e "BiocManager::install('MultiDataSet')"
RUN apt-get update && apt-get install -y \
 
    libbz2-dev \
    liblzma-dev

RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('CTDquerier')"

# select port
 
EXPOSE 80
 
# allow permission for user ‘shiny’ to run
 
RUN sudo chown -R shiny:shiny /srv/shiny-server
  
# install linux programs to enable conversion of ms dos file to unix file
 
# RUN apt-get update && apt-get install -y dos2unix
 
# Change access permissions to shiny-server.sh - did not need this for my purposes
 
RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"] 
 
# run app
 
CMD ["/usr/bin/shiny-server.sh"]

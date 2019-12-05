FROM rocker/r-ver:devel

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    wget

# Download and install shiny server
RUN wget --no-verbose http://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "http://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'shinyjs', 'shinydashboard', 'shinyWidgets', 'shinyBS', 'shinycssloaders'), repos='http://cran.rstudio.com/')" 

# Install R packages that are required
# TODO: add further package if you need!
RUN R -e "install.packages(c('DT', 'ggplot2', 'ggrepel', 'ggseqlogo', 'ggsci', 'pheatmap', 'dplyr', 'stringr', 'knitr', 'reshape2'), repos='http://cran.rstudio.com/')" 

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY app /srv/shiny-server/apps

# Make the ShinyApp available at port 3838
EXPOSE 3838

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh
CMD ["/usr/bin/shiny-server.sh"]
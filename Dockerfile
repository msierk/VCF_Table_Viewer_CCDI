FROM rocker/r-ver:4.5.1

RUN apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc zlib1g-dev curl xtail wget libxt-dev openssl1.1 aptitude && rm -rf /var/lib/apt/lists/* 

RUN echo "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/all/__linux__/focal/latest'), download.file.method = 'libcurl', Ncpus = 4)" >> /usr/local/lib/R/etc/Rprofile.site 

RUN touch init-shiny

RUN mkdir /vcftableviewer 

ADD renv / vcftableviewer

ADD renv.lock /vcftableviewer

WORKDIR /vcftableviewer

RUN R -e 'install.packages("renv"); renv::consent(provided=TRUE)' 

ENV RENV_CONFIG_REPOS_OVERRIDE=https://packagemanager.rstudio.com/all/__linux__/focal/latest 

RUN Rscript -e 'renv::restore()' 

ADD . /vcftableviewer

EXPOSE 3838 

RUN echo "options('shiny.port'=3838,shiny.host='0.0.0.0');OmicCircosShiny::run_app()" >> /init-shiny 

CMD R -e "options('shiny.port'=3838,shiny.host='0.0.0.0');OmicCircosShiny::run_app()" 
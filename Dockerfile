FROM rocker/shiny:4.6.1
ARG RENV_PATHS_CACHE=/root/.cache/R/renv
ENV "RENV_PATHS_CACHE"="${RENV_PATHS_CACHE}"
RUN apt-get update -y && apt-get install -y  cmake make libuv1-dev pandoc libcurl4-openssl-dev libssl-dev libbz2-dev liblzma-dev libxml2-dev zlib1g-dev libicu-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(renv.config.pak.enabled = FALSE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("renv", version = "1.1.5")'
WORKDIR /srv/shiny-server/
COPY . /srv/shiny-server/
RUN R -e 'renv::restore()'
EXPOSE 3838
CMD R -e 'shiny::runApp("/srv/shiny-server",host="0.0.0.0",port=3838)'

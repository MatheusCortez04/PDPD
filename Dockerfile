
FROM rocker/tidyverse:4.5.2


USER root


RUN apt-get update && apt-get install -y libglpk-dev && \
    rm -rf /var/lib/apt/lists/*


RUN R -e "install.packages(c('here','BiocManager','reportROC'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('diffuStats', 'igraphdata', 'igraph'), ask=FALSE)"

WORKDIR /home/rstudio

COPY .  .
RUN chown -R rstudio:rstudio /home/rstudio
USER rstudio
CMD ["R"]
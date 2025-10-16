
FROM rocker/rstudio:latest


USER root


RUN apt-get update && apt-get install -y libglpk-dev && \
    rm -rf /var/lib/apt/lists/*


RUN R -e "install.packages(c('here', 'Matrix', 'BiocManager'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('diffuStats', 'igraphdata', 'igraph'), ask=FALSE)"

WORKDIR /home/rstudio

COPY ./src .
CMD ["Rscript", "PDPD.R"]
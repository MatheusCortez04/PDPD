
FROM bioconductor/bioconductor_docker:devel-r-4.6.0

USER root


RUN apt-get update && apt-get install -y libglpk-dev && \
    rm -rf /var/lib/apt/lists/*


RUN R -e 'install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))'
RUN R -e "pak::pkg_install(c('tidyverse','here', 'reportROC', 'httr2', 'diffuStats', 'igraphdata', 'igraph', 'clusterProfiler', 'org.Hs.eg.db'))"
WORKDIR /home/rstudio

COPY .  .
RUN chown -R rstudio:rstudio /home/rstudio
USER rstudio
CMD ["R"]
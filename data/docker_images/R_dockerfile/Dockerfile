FROM rocker/tidyverse:4.0.5

ARG WHEN

ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /project
COPY renv.lock renv.lock
RUN R -e 'renv::restore()'

RUN R -e "install.packages('e1071', dependencies=TRUE)"


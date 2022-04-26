FROM rocker/tidyverse:4.0.5

ARG WHEN

ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /project

RUN R -e "install.packages('C50', dependencies=TRUE)"
RUN R -e "install.packages('caret', dependencies=TRUE)"

FROM rocker/tidyverse:4.0.5

ARG WHEN

RUN R -e "install.packages('radomForest',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('caret',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('doParallel',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# TODO CHANGE
#RUN mkdir /home/analysis
#COPY myscript.R /home/analysis/myscript.R

# docker build -t staphb/piggy:1.5 .
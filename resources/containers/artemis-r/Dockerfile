FROM rocker/tidyverse:latest

# Metadata
LABEL author="Okawa"
LABEL version="1.0"
LABEL description="Custom R Tidyverse image with specific versions of argparse and data.table."

# Environment variables
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN R -e "devtools::install_version('argparse', version = '2.3.1', repos = 'http://cran.us.r-project.org')" \
    && R -e "devtools::install_version('data.table', version = '1.17.8', repos = 'http://cran.us.r-project.org')" \
    && rm -rf /tmp/downloaded_packages/

CMD ["R"]
#docker build -t artemis-r:1.0.0 .

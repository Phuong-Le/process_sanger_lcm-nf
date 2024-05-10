FROM ubuntu:24.04
LABEL maintainer="Phuong Le" 
LABEL maintainer.email="al35@sanger.ac.uk"



RUN apt-get update \
	&& apt-get install -yq build-essential wget default-jdk \
	&& wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda.sh \
	&& /bin/bash ./miniconda.sh -b -p /opt/conda \
	&& rm -f ./miniconda.sh \
	&& /opt/conda/bin/conda init bash \
    && apt-get -y install git-all
ENV PATH=/opt/conda/bin:$PATH

RUN conda install -y -c conda-forge \
    r-base=4.3.3 \
    r-devtools \ 
    r-biocmanager \
    r-data.table \
    r-vgam \
    r-ape \
    r-ggplot2 
RUN conda install -y -c bioconda \
    bioconductor-ggtree 
RUN R -e "devtools::install_git('https://github.com/NickWilliamsSanger/treemut')"


# testing that everything is installed correctly
RUN R -e "library(data.table); \
    library(VGAM); \
    library(ape); \
    library(ggplot2); \
    library(ggtree); \
    library(treemut); \
    print('all packages are ready to be used')"

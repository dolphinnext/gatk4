FROM ubuntu:xenial

RUN apt -y update

RUN apt -y install \
	git \
	wget \
	autoconf \
	automake \
	make \
	gcc \
	perl \
	zlib1g-dev \
	bzip2 \
	libbz2-dev \
	xz-utils \
    liblzma-dev \
	curl \
    libcurl4-openssl-dev \
	libssl-dev \
	ncurses-dev \
	graphviz \
    unzip \
    zip

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

COPY environment.yml /
RUN . /opt/conda/etc/profile.d/conda.sh && \ 
    conda activate base && \
    conda update conda && \
    conda install -c conda-forge mamba && \
    mamba env update --file /environment.yml --prune && \
    mamba clean -a

RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/bin:/opt/conda/condabin:$PATH
ENV CONDA_EXE /opt/conda/bin/conda
ENV JAVA_HOME /opt/conda
ENV CONDA_PREFIX /opt/conda
ENV JAVA_LD_LIBRARY_PATH /opt/conda/lib/server
ENV CONDA_SHLVL 1
ENV CONDA_PROMPT_MODIFIER (base)
ENV CONDA_PYTHON_EXE /opt/conda/bin/python
ENV CONDA_DEFAULT_ENV base

# R Packages Installation
COPY install_packages.R /
RUN Rscript /install_packages.R


#############################
# CentOS-based Docker build #
#############################

## The Base OS and APBS build are taken from the APBS Centos Dockerfile

###################
# Install Base OS #
###################

FROM centos:7 AS apbs_base_centos

RUN yum -y install \
        epel-release \
        centos-release-scl && \
    yum-config-manager --enable rhel-server-rhscl-7-rpms && \
    yum -y install devtoolset-11 && \
    source scl_source enable devtoolset-11 && \
    yum -y install \
        unzip \
        cmake3 \
        git \
        arpack-devel \
        arpack-static \
        f2c \
        eigen3-devel \
        boost-devel \
        python3-pip \
        python3-devel \
        openblas-serial64 \
        openblas-static \
        openblas-devel \
        lapack64-devel \
        lapack64-static \
        suitesparse-devel \
        suitesparse-static \
        && \
    ln -s /usr/bin/cmake3 /usr/local/bin/cmake

##########################
# Install METIS for APBS #
##########################

FROM apbs_base_centos

RUN echo "${HOME}"

RUN source scl_source enable devtoolset-11 && \
    git clone https://github.com/Electrostatics/metis-fork.git metis-5.1.0 && \
    cd metis-5.1.0 && \
    make config prefix=/usr/local && \
    make install

################
# Install APBS #
################

FROM apbs_base_centos

ADD /basechem/common/apbs /tmp_source

ARG BLA_VENDOR="OpenBLAS"
ARG INSTALL_DIR=/usr/local
ARG RELEASE_TYPE=Debug
ARG ENABLE_BEM=ON
ARG FETK_VERSION=57195e55351e04ce6ee0ef56a143c996a9aee7e2
ARG PYTHON_VERSION="3.6"
ARG PYTHON3_INCLUDE_DIR="/usr/include/python3.6m"
ARG PYTHON3_LIBRARY="/usr/lib64/libpython3.6m.so"
ARG MAKEJOBS="-j"

RUN source scl_source enable devtoolset-11 && \
    cd /tmp_source && \
    mkdir build && cd build && \
    cmake \
      -DCMAKE_INSTALL_INCLUDEDIR="include" \
      -DBUILD_DOC=ON \
      -DAPBS_STATIC_BUILD=ON  \
      -DBUILD_TOOLS=ON \
      -DCMAKE_BUILD_TYPE=$RELEASE_TYPE \
      -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
      -DENABLE_PYGBE=ON \
      -DENABLE_BEM=ON \
      -DENABLE_iAPBS=ON \
      -DENABLE_GEOFLOW=ON \
      -DENABLE_OPENMP=OFF \
      -DENABLE_PBAM=OFF \
      -DENABLE_PBSAM=OFF \
      -DENABLE_PYTHON=OFF \
      -DENABLE_TESTS=ON \
      -DFETK_VERSION=${FETK_VERSION} \
      -DGET_NanoShaper=ON \
      -DPYTHON_VERSION="${PYTHON_VERSION}" \
      -DPython3_INCLUDE_DIR=${PYTHON3_INCLUDE_DIR} \
      -DPython3_LIBRARY=${PYTHON3_LIBRARY} \
      .. && \
    make ${MAKEJOBS} install


######################################
# Set up base linux env for Basechem #
######################################
    
ARG GITHUB_TOKEN

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN yum install wget -y 

RUN mkdir -p ~/miniconda3 && \ 
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --no-check-certificate -O ~/miniconda3/miniconda.sh && \ 
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 && \ 
    rm -rf ~/miniconda3/miniconda.sh && \ 
    ~/miniconda3/bin/conda init bash

RUN source ~/.bashrc && conda config --add channels conda-forge --add channels bioconda
RUN source ~/.bashrc && conda config --remove channels defaults

RUN ln -sf /usr/share/zoneinfo/US/Pacific /etc/localtime

RUN yum update -y && yum install -y \
    unzip \
    vim-enhanced \
    python-devel \
    epel-release \
    bzip2

# Install git >= 2.18 for use in GitHub workflows
RUN yum -y install https://packages.endpointdev.com/rhel/7/os/x86_64/endpoint-repo.x86_64.rpm \
    && yum install -y git

RUN yum install -y \
    libXcomposite libXcursor libXi libXtst libXrandr alsa-lib mesa-libEGL libXdamage mesa-libGL libXScrnSaver libuv \
    libxml2-devel xmlsec1-devel xmlsec1-openssl-devel libtool-ltdl-devel xmlsec1 gcc \
    python-psycopg2 nc \
    supervisor sudo

# Set up a conda environment. Use conda-libmamba-solver as the dependency solver to speed up installation of dependencies.
RUN source ~/.bashrc && conda update -n base conda -y \
    && conda install -n base conda-libmamba-solver \
    && conda config --set solver libmamba \
    && conda init bash

# lxml is pinned to 4.9.3 because 5.0.0 is incompatable with python3-saml (at least through v1.16.0), which is the package we use for Okta SAML integration.
# This issue is described here: https://github.com/SAML-Toolkits/python3-saml/issues/389
RUN source ~/.bashrc && conda create --name cenv python=3.10 scikit-bio=0.5.9 rdkit=2022.09.4 lxml==4.9.3 nodejs

RUN echo "conda activate cenv" >> ~/.bashrc
ENV VIRTUAL_ENV /usr/local/envs/cenv
ENV PATH /usr/local/envs/cenv/bin:$PATH

# Make user home
RUN useradd -U -m -s /bin/bash -G wheel mniadmin

# Install dtxwrapper for Dotmatics access, pinned commit to reduce possibility of breaking
# ** Remove this line if building locally and not part of Denali ** 
RUN curl -L https://github.com/denalitherapeutics/dtxwrapper/archive/3e059d1bca241fc1c268429544872c27a56b1155.tar.gz -H "Authorization: token $GITHUB_TOKEN" --output /tmp/dtxwrapper.tar.gz && \
    cd /tmp && tar -zxvf dtxwrapper.tar.gz && rm -rf dtxwrapper.tar.gz && mv dtxwrapper-3e059d1bca241fc1c268429544872c27a56b1155 dtxwrapper

RUN /root/miniconda3/envs/cenv/bin/python -m pip install -e /tmp/dtxwrapper -vv


####################################
### Set up Basechem dependencies ###
####################################

# Install pango; needed for weasyprint
RUN yum install -y pango

# Install rDock + dependencies
RUN yum group install "Development Tools" -y
RUN yum install cppunit cppunit-devel popt-devel wget -y
RUN wget -P /opt/ https://github.com/CBDD/rDock/releases/download/v0.22/rDock_2022_v0.tar.gz && cd /opt && tar -xvzf rDock_2022_v0.tar.gz && rm rDock_2022_v0.tar.gz
RUN cd /opt/rDock_2022_src/build/ && make clean && make linux-g++-64
ENV RBT_ROOT=/opt/rDock_2022_src
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RBT_ROOT/lib
ENV PATH=$PATH:$RBT_ROOT/bin

# Install Mayachemtools
RUN wget -P /opt/ http://www.mayachemtools.org/download/mayachemtools.tar.gz && cd /opt && tar -zxvf mayachemtools.tar.gz && rm mayachemtools.tar.gz
# Install Mayachemtools dependencies
RUN source ~/.bashrc && conda install -n cenv seaborn psi4 dftd3-python gcp-correction -c conda-forge

# Install LS Align
RUN wget -P /opt/LS-align https://seq2fun.dcmb.med.umich.edu//LS-align/LSalign_src.tar.gz && cd /opt/LS-align && tar -xvzf LSalign_src.tar.gz && cd src  && g++ LSalign.cpp -o LSalign && rm ../LSalign_src.tar.gz
RUN  source ~/.bashrc && conda install --name cenv -c openbabel openbabel

# Install ESP-DNN environment; Python 2.7
COPY ./basechem/common/ESP_DNN/environment.yml . 
RUN  source ~/.bashrc && conda env create --name esp-dnn-env -f environment.yml

# Install Toklat 
ADD ./basechem/common/toklat /opt/toklat
RUN source ~/.bashrc && cd /opt/toklat && pip install -e .

# Download MMPDB
RUN source ~/.bashrc && wget -P /opt/ https://github.com/rdkit/mmpdb/archive/refs/tags/v3.1.tar.gz && cd /opt && tar -zxvf v3.1.tar.gz && rm v3.1.tar.gz

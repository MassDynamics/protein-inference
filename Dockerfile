FROM ubuntu:18.04

USER root
ENV DEBIAN_FRONTEND noninteractive 

RUN echo "Installing non python dependencies"
RUN apt-get update
RUN apt-get install -y unzip gcc build-essential
RUN apt-get install -y cmake g++ autoconf qtdeclarative5-dev patch libtool make git qtbase5-dev libqt5svg5-dev libqt5opengl5-dev automake libqtwebkit-dev
RUN apt-get install -y libboost-regex-dev libboost-iostreams-dev libboost-date-time-dev libboost-math-dev libsvm-dev libglpk-dev libzip-dev zlib1g-dev libxerces-c-dev libbz2-dev
RUN apt-get install -y wget

RUN echo "Installing Anaconda3"
RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -P /tmp
RUN chmod a+w /opt
RUN bash /tmp/Anaconda3-2020.11-Linux-x86_64.sh -b -p /opt/anaconda3
RUN ln -s /opt/anaconda3/bin/* /usr/local/bin
RUN conda install -y -c anaconda cython

RUN echo "Installing AWS CLI..."
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

RUN echo "Setting Up conda environment"
RUN conda init bash
COPY ./requirements.txt .
RUN conda config --append channels conda-forge
RUN conda config --append channels scikit-learn
RUN conda create --name protein-inference python=3.8 --file requirements.txt

FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY ./requirements.txt .
RUN conda config --append channels conda-forge
RUN conda config --append channels scikit-learn
RUN conda create -y --name protein-inference python=3.8 --file requirements.txt

# Make RUN commands use the new environment:
RUN echo "conda activate protein-inference" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Demonstrate the environment is activated:
RUN echo "Make sure pyvis is installed:"
RUN python -c "import pyvis"

# install protein inference package
RUN echo "Installing Protein Inference Python Package"
COPY ./setup.py .
ADD protein_inference/ ./protein_inference
RUN python setup.py install



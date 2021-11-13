#!/bin/bash

set -e -x

conda env remove --name protein-inference
conda env create --name protein-inference --file requirements.txt python=3.8.5 

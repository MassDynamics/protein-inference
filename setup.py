#!/usr/bin/env python

import setuptools
from distutils.core import setup
from setuptools import find_packages

setup(name='protein_inference',
      version='0.2.4',
      description='Protein Inference Library for Network based Inference',
      long_description=
      "This codebase is being developed at Mass Dynamics https://www.massdynamics.com/ to facilitate work on the [Protein Inference Problem](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-S16-S4). More here: https://github.com/MassDynamics/protein-inference",
      author='Joseph Bloom',
      author_email='joseph@massdynamics.com',
      url='https://www.massdynamics.com/',
      packages= find_packages(exclude=['example_code','example_data']),
      install_requires=[
          "pandas>=1.1.4",
          "networkx>=2.3",
          "pyvis>=0.1.8.2",
          "matplotlib>=3.2.2",
          "upsetplot>=0.4.1",
          "seaborn>=0.11.0",
          "plotly_express>=0.4.1",
          "plotly>=5.7.0",
          "streamlit>=1.8.1",
          "streamlit-plotly-events>=0.0.6")
#!/usr/bin/env python

import setuptools
from distutils.core import setup
from setuptools import find_packages

setup(name='protein_inference',
      version='0.2.1',
      description='Protein Inference Library for Network Based Inference',
      author='Joseph Bloom',
      author_email='joseph@massdynamics.com',
      url='https://www.massdynamics.com/',
      packages= find_packages(),
      install_requires=[
          "pandas>=1.1.4",
          "networkx>=2.3",
          "pyvis>=0.1.8.2",
          "matplotlib>=3.2.2",
          "upsetplot>=0.4.1",
          "seaborn>=0.11.0",
          "plotly_express>=0.4.1",
          "plotly>=4.10.0"])
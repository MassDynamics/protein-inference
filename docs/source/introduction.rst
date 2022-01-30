============
Introduction
============

The goal of this project is to build a graph based approach to performing protein inference calculations, 
that makes inferences about which proteins likely emitted which peptides (that were observed via PSM's),
thus creating an interpretable basis for protein quantification.

The protein inference problem is the challenge of inferring the existence of proteins in a sample based 
on peptide spectral matches that result from search algorithms such as percolator. 

This project makes extensive use of Networkx and pandas. It is written in the hope that further work can be built 
using this framework.
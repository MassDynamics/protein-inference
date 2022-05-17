
============
Tutorial
============

The protein_inference library is meant to facilitate scientific thinking
around protein inference including algorithms that perform protein grouping, 
peptide assignment and ultimate protein scoring/
inference. It is recommended to understand this problem well before
proceeding to make use of the package by reading this `paper`_ (as a start.)

.. _paper: https://academic.oup.com/bib/article/13/5/586/415393


Using the Package as is:
------------------------

This package can be used with default settings via the commandline or 
via the ProteinInferenceRunner() class.

Bash

>>> python -m <virtualenv_name>/lib/<python_ver>/site-packages/protein_inference.main \
    --output-directory path/to/out/folder \
    --target-path path/to/your/target.psms.txt \
    --decoy-path path/to/your/decoy.psms.txt

Python 

.. code-block:: python

    import protein_inference as pi

    target_path = "path/to/your/target.psms.txt"
    decoy_path =  "path/to/your/decoy.psms.txt"
    output_directory = "path/to/your/output_directory"
    pi.ProteinInferenceRunner().run(target_path, 
                                decoy_path,
                                output_directory)

Some generic solvers are available in inference.scorers which can 
be passed to the ProteinInferenceRunner().run() method via the scoring method
argument.

For example:

.. code-block:: python

    import protein_inference as pi

    target_path = "path/to/your/target.psms.txt"
    decoy_path =  "path/to/your/decoy.psms.txt"
    output_directory = "path/to/your/output_directory"
    pi.ProteinInferenceRunner().run(target_path, 
                                decoy_path,
                                output_directory,
                                scoring_method = pi.inference.scorers.PEPProductScorer())



Using the package to visualize PSM Networks:
--------------------------------------------

One of the most exciting and enjoyable components of this package is the network visualizations
of psm networks. By default, all target problem networks (psm networks) are pickled inside a list
when the main workflow is run. 

To access and visualize these networks, locate the `target_networks.p` pickle file and load 
it into an ipython notebook such as jupyter. 

Then you can find molecules and visualize networks by running code like the code below. The "by" argument
to Network grapher allows you to visualize networks by allocated peptide and protein "status",
by allocated "groups" or by "score".


Here I used the default status visualization:

.. code-block:: python

    import protein_inference as pi
    import pickle

    target_networks = pickle.load(open("example_data/outputs/target_networks.p", "rb"))
    large_target_networks = [pn for pn in target_networks if len(pn.get_proteins())>1]
    NetworkGrapher().draw(large_target_networks[0], by="status", 
                    name="example_status", size = [400,600])
    
.. raw:: html
	:file: example_status.html


Using the package to implement a Protein Inference Strategy:
------------------------------------------------------------

If you are developing or would like to try implementation of a protein inference algorithm, 
it is highly recommended that you understand the protein inference workflow as it is currently
implemented. This can be easily seen in protein_inference.ProteinInferenceRunner().get_output()
which shows the steps taken in parallel for each of the target and decoy datasets. 

These steps are (classes in brackets):

* Table Preprocessing (PSMsPreprocessor)
* Conversion of the PSM table to one large PSM Network (PSMsNetworkGenerator)
* Splitting the large PSM Network (PSMsNetworkSplitter)
* Tagging Unique peptides and unique evidenced Proteins (UniquenessTagger)
* Solving Networks (ie: Scoring and/or adding addition tags) (see in solvers.py or reprisal.GreedyAlgorithm)
* Merging Proteins with identical PSM connections (ProteinMerger)
* Generating Output Tables (TableMaker)

After scores have been calculated for target and decoy tables, ProteinInferenceRunner will
call FalseDiscoveryRateCalculator to estimate FDR's or q-values. 

There are two ways to develop a protein inference strategy using protein_inference: 


#. Use the REPRISAL workflow directly, substituting a solver such as those in solvers.py for the protein scoring calculations. 
#. Write a bespoke ProteinInferenceRunner() using any subset/combination of the functions inside the package. 

Running Unit Tests
------------------

While developing, the core functionality of the workflow is protected by unit tests. Unit tests are
automatic tests developed during with the main functionality that check it works before and after 
edits are made. 

To run the unit tests from the root directory.

Bash

>>> ./scripts/acceptance.sh

In order to maintain the ongoing integrity and extemsibility of the code, it is an expectation of contributors that
merge request code can pass all of the previous build tests and new functionality be protected with new build tests.

To write build tests simply follow the patterns from existing code. 
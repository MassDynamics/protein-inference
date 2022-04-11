====================
Preparing Input Data
====================

In order to perform protein inferece in a bottum up proteomics LC-MS/MS pipeline, we must  first 
perform the following steps:

1. Conversion of .RAW files (or other vendor-specific Mass Spectrometer output files) to mgf or other search-engine compatible formats.

2. Performing an MS/MS search using tools such as Comet or Mascot. These tools read MS/MS fragmentation spectra return a list of matches between these spectra and theoretical spectra produced from a hypothesis fasta file (sequence database). Each match between an MS/MS spectrum and a possible peptide is called a PSM (peptide-spectrum match).

3. Scoring can then be performed to more carefully distinguish between good and poor matches. The most common method is called the "target-decoy approach" which involves searching a "decoy" database of peptides known not to be present in the sample. 

4. Finally, once PSM's have been scored, this information can be aggregated to the protein level, answering the question "Which proteins are present in my sample?".

Below, I have provided an example pipeline (designed for identification of peptides in Thermo-Fischer .Raw files) using docker containers to remove as much complexity as possible. 

The example data is the iPRG2016 Protein Inference Benchmark dataset.

The software I will use for each stage is:

1. msconvert
2. comet via crux
3. percolator via crux
4. protein inference (this python package)

Converting Raw files with msConvert
-----------------------------------

First, you will need to convert your raw files to mzML files which can be used as an input for each. 

.. code-block:: bash

    # Get the msconvert Docker Container
    docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

    experimentfolder=${1} # will become the volume in the docker command
    relative_upload=${2} # relative to the experiment folder, where are the mzML's

    docker run -it --rm -e WINEDEBUG=-all -v $experimentfolder:/data \
        chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/${2}/*.raw \
        --filter "peakPicking true" -o /data/ms_convert_output



Performing a search with Crux Comet
-----------------------------------

Before using crux, you will need to build a docker container with the latest version of crux.

.. code-block:: bash

    # installing crux
    git clone https://github.com/crux-toolkit/crux-toolkit.git
    # update the docker file with credentials

    cd crux-toolkit

    # you will need to add your github information to the docker file first!
    docker-compose build # build it 
    docker-compose up # test it 

    # get help 
    docker run -it --rm  \
        cruxtoolkit_crux crux comet /data/


When you have successfully installed a crux docker image, you can then use it to perform a search with compatible. 

Notes:

- Chosen parameters are not intended as recommendations.
- You may wish to pass multiple mz_Ml files to comet simultaneously. 

.. code-block:: bash

    experimentfolder=${1} # will become the volume in the docker command
    mz_ML_relative_location=${2} # relative to the experiment folder, where are the mzML's
    fasta_relative_location=${3} # relative to the experiment folder, where is the fasta database

    docker run -it --rm -v $experimentfolder:/data \
        cruxtoolkit_crux crux comet \
        $mz_ML_relative_location \
        $fasta_relative_location \
        --decoy_search 1 \
        --output_percolatorfile 1 \
        --output-dir /data/comet \
        --peptide_mass_tolerance 20 \
        --peptide_mass_units 2 \
        --isotope_error 0 \
        --fragment_bin_tol 0.01 \
        --allowed_missed_cleavage 1 \
        --peptide_length_range "5 63" \
        --overwrite T


PSM Scoring with Percolator
---------------------------

Now you can take the output pep.xml files (which also contain the decoy search results),
and use percolator for PSM scoring.

Notes:

- Chosen parameters are not intended as recommendations.
- You may wish to pass multiple mz_Ml files to comet simultaneously. 


.. code-block:: bash

    experimentfolder=${1} # will become the volume in the docker command
    pep_xml_relative_location = ${2} # relative to the experiment folder, where is comet output.

    # run percolator on 1 file ... 
    docker run -it --rm -v $experimentfolder:/data \
        cruxtoolkit_crux crux percolator \
        $pep_xml_relative_location
        --verbosity 30 \
        --output-dir /data/percolator \
        --enzyme trypsin \
        --only-psms true \
        --overwrite T


Wrapping Up
-----------

Congratulations! If you've gotten this far, you have:

1. Built docker images for msconvert and crux 
2. Converted .raw files mzMl search-ready files. 
3. Performed search using Comet. 
4. Perfomed PSM scoring with Percolator. 

Most importantly, you are now ready to perform Protein Inference with your chosen method!

This python package contains a framework for protein inferencei in python, including
a novel algorithm "RePrISAl" (Reprisal) which performs recursive assignment of peptides to 
proteins based on uniqueness of evidence and score to enable interpretable protein inference. 

If you have found this tutorial useful or had any trouble, or would just like to chat, please let me know @ joseph@massdynamics.com!


Identification Pipeline Resources:
----------------------------------

* msConvert
* Crux 
* Comet 
* Percolator

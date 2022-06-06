# Protein Inference - Mass Dynamics (MD)

This codebase is being developed at Mass Dynamics https://www.massdynamics.com/ to facilitate work on the [Protein Inference Problem](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-S16-S4). 

It has been developed using  the output of [crux comet/percolator](https://crux.ms/commands/percolator.html) (eg: percolator.target.psms.txt), but requires only that psm tables contain q-value and score columns, although more information about psms could theoretically be integrated into scoring methods. 

This codebase could be used for the following purposes:

    * to develop or benchmark protein inference algorithms using shared data processing and  parallelization infrastructure
    * to visualize the peptide-protein networks that are being reasoned about by a protein inference algorithm
    * to opensource an algorithm that is implimented inside this package

The REPRISAL algorithm, currently unpublished, was the first novel algorithm implemented inside the protein_inference framework and is contained in the distinct subpackage "reprisal". 

REPRISAL stands for "REcursive PRotein Inference Scoring ALgorithm" and works by iteratively scoring proteins by their associated peptides and grouping all associated and unallocated peptides with the highest scoring protein at each iteration. FDR's are then calculated as per the score distributions for proteins in the target and decoy tables. 

## Installation 

```
pip install protein-inference
```

## Build Docs

```
sphinx-build -b html docs/source/ docs/build
```

## Development Set up

```
conda create -y --name protein-inference python=3.8 --file requirements.txt
```

## Dependencies

Many of these are just used in the benchmarking suite for visualization. 

[pandas](https://github.com/pandas-dev/pandas)

[networkx](https://networkx.org/)

[pyvis](https://pyvis.readthedocs.io/en/latest/index.html)

[matplotlib](https://matplotlib.org/)

[plotly](https://plotly.com/python/)

[seaborn](https://seaborn.pydata.org/)

[upsetplot](https://pypi.org/project/UpSetPlot/)

## Quick Start

If you have a psms list from percolator (or otherwise) ready to go, you can run the default pipeline easily:

```python
import protein_inference as pi

target_path = "path/to/your/target.psms.txt"
decoy_path =  "path/to/your/decoy.psms.txt"
output_directory = "path/to/your/output_directory"
pi.ProteinInferenceRunner().run(target_path, 
                            decoy_path,
                            output_directory)

```

You can also run the package via the commandline if you can specify the path to the python directory where you have protein_inference installed. 

```bash
python -m protein_inference.main \
  --output-directory path/to/out/folder \
  --target-path path/to/your/target.psms.txt \
  --decoy-path path/to/your/decoy.psms.txt

```

## Streamlit App

If you'd like to play around with some processed data/visualization, install streamlit and run the protein inference playground app.

```bash
pip install streamlit
cd streamlit
streamlit run pi_playground.py -- --folder /path/to/your/pi/output/
```
You can also start the app without a path to an upload folder, and upload your output files after the fact.

## Using Docker:

A docker file has been created which creates an environment to run the protein_inference
package via a command line interface. 

You can build the image from the docker file as so: 

```bash
docker build --pull --rm -f "Dockerfile" -t proteininference:latest "." 
```

Then you can run the main file from the command line by using a /data volume and passing arguments with paths to files in that folder.

```bash
docker run -it --rm -v /data:/data \
    proteininference conda run -n protein-inference -c python -m protein_inference.main 
    --output-directory /data/output \
    --target-path /data/input/percolator.target.psms.txt \
    --decoy-path /data/input/percolator.decoy.psms.txt
```

## Future Plans: 

The benchmarking suite is currently undocumented as it is still being developed. I'll get around to this as I work towards publishing REPRISAL as a paper. 

While not currently in development, possible additions to this codebase are:
  - Writing file converters to accept mzIdentML
  - Integrating Baysian methods such as Fido or EPIFANY.
  - Integrating an optimisation method, maybe using a package like PULP. 
  - Further development of benchmarking suite to allow easy comparison of existing algorithms. 

## ASMS 2022 Poster

![ASMS 2022 Poster](ASMS2022_Poster.png)

Feel free to download the png or pdf of the ASMS 2022 poster on protein inference in python. 

<a href="ASMS2022_Poster.pdf.js" download>Click to Download the poster pdf</a>

NB: Table 1 was removed when trying to improve the app asthetics. I've provided it here for reference


![Table 1](table1.png)

## App Link!

https://share.streamlit.io/massdynamics/protein-inference/main/pi_playground.py
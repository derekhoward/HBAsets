HBASets
=======

HBASets is a project for creating molecular maps of the brain starting with a query of a gene set of interest.

## What is this tool?
Brain diseases are often due to variations in multiple genes (polygenic disorders) and their interactions with environmental factors.

There is a growing amount of open data describing gene expression in the brain. For example, the [Allen Brain Atlas](http://www.brain-map.org/) provides a number of highly comprehensive brain atlases in the human, mouse and monkey. However, these valuable resources are underused for the analysis of polygenic brain disorders because the data is not easily accessible beyond the level of a single gene.

This tool aims to make use of the substantial open data of gene expression in the brain to facilitate accessible, rapid and custom data mining of open brain transcriptome data (across time, anatomy, species and celltypes).

## Who is this for?
Many scientific studies produce lists of genes that are differentially expressed in a brain disease or where genetic variants are associated with a mental illness. Often, these gene sets are derived from unbiased genome-wide studies and may not have been previously characterized.

This tool will further analyses and help answer questions such as:
- Are these genes most expressed in childhood?
- Are they expressed in a specific brain area?
- Are they expressed in neurons or glia?
- Do the above answers agree in mouse, monkey and human brain?


## Data
To start with, you can download the [normalized microarray datasets](http://human.brain-map.org/static/download) of gene expression from 6 adult human brains that was released by the Allen Brain Atlas. Place the raw data in /HBAsets/data/raw/ and you can run the make_dataset.py preprocessing script to generate the processed expression matrix which will be used for analysis.

Alternatively, if you would like to go straight to the analysis steps download the preprocessed data from GDrive [here](https://drive.google.com/drive/folders/0BzB3t6aSc9bDdDc2MEU3MjROY0E?usp=sharing) and place the files directly in /HBAsets/data/processed


## How can I get involved?
HBAsets is openly developed and welcomes contributors. Check out the [Roadmap](https://github.com/derekhoward/HBAsets/wiki) and the [contributing guidelines](CONTRIBUTING.MD) to get an idea of how to start. Your thoughts, questions, bugs, suggestions, etc are valuable to us so feel free to take part in discussion in the [issues](https://github.com/derekhoward/HBAsets/issues).

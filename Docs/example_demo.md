# Worked example

This guide will walk through the main Jupyter-Notebook of the software - lrRNAseq_GAST_visualiser.ipynb, cell-by-cell, explaining in detail, how to use each cell, what is being done, and the data structures generated in the process. 

For best experience, Linux-based operating systems are advised; this guide will assume you have a Linux-based operating system. This is mostly important when using BLAST.

## Step 1 - Download lrRNAseq-GAST

Download the entire lrRNAseq-GAST repository. This can be done by clicking on the green "<> Code" drop down menu button on the main page of lrRNAseq-GAST repository, then clicking "Download ZIP", and extracting the ZIP into the directory of choice.

An alternative method is to use git:
```
git clone https://github.com/Maxim-Karpov/lrRNAseq-GAST.git
```

## Step 2 - Install Anaconda Python distribution

For novice users it is recommended to install Anaconda distribution platform which contains Python, Jupyter-Notebook, and some required Python libraries such as Pandas and Matplotlib, and adds Python to your shell environment. Follow the instructions on the official Anaconda website - https://www.anaconda.com/download. 

Once installed, open Jupyter-Notebook through Anaconda Navigator and launch lrRNAseq-GAST by opening lrRNAseq_GAST_visualiser.ipynb. The full list of dependencies includes: Pandas, Matplotlib, Numpy, Regex, Biopython, OS, ast, operator. With Anaconda, you should only need to manually install biopython and seqtk:
```
conda install anaconda::biopython
conda install bioconda::seqtk
```

## Step 3 - Check all imports

To check whether all dependencies are present in your Python environment, and the supplementary functions can be imported successfully, run the first cell of the notebook. The bioinformatic_functions.ipynb must be in the same directory as lrRNAseq_GAST_visualiser.ipynb for this to work. In case some dependencies are missing, consult the internet on their installation.

```
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import matplotlib.patches as patches
import matplotlib
from Bio import SeqIO
from matplotlib.ticker import MaxNLocator
from ast import literal_eval
from operator import itemgetter
import warnings
warnings.filterwarnings("ignore")
%run "./bioinformatic_functions.ipynb"
```

# Step 4 - Align long read RNAseq transcripts to a genomic region

Install BLAST in your shell environment:

```
conda install bioconda::blast
```

Make BLAST database of your genomic sequence. Here we're using Ef318A10.fst (provided in the repository), which is a Bacterial Artificial Chromosome sequence of a genomic locus containing an earthworm metallothionein gene, about 100 kb in length. In your shell environment, whilst in the directory containing Ef318A10.fst genomic sequence, run:

```
makeblastdb -in Ef318A10.fst -dbtype nucl -parse_seqids -out 318A10_db -title 318A10_db
```

Align long read RNAseq transcripts to a genomic region using BLAST with alignment parameters of your choice. Example transcript file "all_transcripts.fasta" was supplied in the lrRNAseq-GAST repository.

```
blastn -query /home/maxim/Projects/BACs/publication/lrRNAseq_GAST/all_transcripts.fasta -db 318A10_db -max_hsps 10000 -task blastn -word_size 4 -evalue 1e-2 -num_threads 10 -out lrRNAseq_to_Ef318A10.out -outfmt '6 qseqid sseqid slen qlen qstart qend sstart send length mismatch gapopen pident evalue bitscore'
```

Important: lrRNAseq-GAST assumes you will use -outfmt '6 qseqid sseqid slen qlen qstart qend sstart send length mismatch gapopen pident evalue bitscore' BLAST output format. 


-word_size 4, -evalue 1e-2 will yield highest sensitivity, and will detect very short exons, but may take unreasonably long to run, and produce many false positives. If such sensitivity is really desired, it is recommended to run BLAST for some managable period of time (e.g. 10 minutes), and interrupt the command; this should still work given your transcripts of interest in your lrRNAseq dataset are not rare / poorly expressed / of low abundance.

-word_size 8, -evalue 1e-5 may be managable to run to the end.

# Step 5 - Import BLAST alignments to your notebook

In the following notebook cell, provide path to your BLAST output. For example:

```
TRANSCRIPT_TO_GENOME_BLAST_PATH = "/home/maxim/Projects/BACs/publication/lrRNAseq_GAST/lrRNAseq_to_Ef318A10.out"
all_blast = read_blast(TRANSCRIPT_TO_GENOME_BLAST_PATH)
```

Then run the cell. This parses the BLAST output to a Pandas dataframe using the read_blast() function. The 

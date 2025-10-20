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

Then run the cell. This parses the BLAST output to a Pandas dataframe using the read_blast() function. The all_blast dataframe contains the BLAST output with some additional columns - sorientation (orientation of subject), q_mid (alignment midpoint on query), s_mid (alignment midpoint on subject).

![Screenshot](./images/all_blast.png)

# Step 6 - Calculate transcript filtration variables

In this step, for each transcript, maximal distance between two alignment regions (maximal intron gap), and transcript coverage (how much of the transcript was covered by all of its alignments to genome), will be calculated for future filtering purposes. Run the cell.

```
grouped_blast = all_blast.groupby("qseqid")

all_blast["Tr_coverage"] = 0
cov_dict = {}
dist_dict = {}
for tr, temp_df in grouped_blast:
    qlen = temp_df["qlen"].max()
    temp_df = temp_df.sort_values(by="qstart").reset_index(drop=True)
    q_overlaps, q_cov = calc_overlaps(temp_df, "q") 
    q_overlaps["Length"] = q_overlaps["Right"] - q_overlaps["Left"] 
    coverage = 100 / qlen * q_overlaps["Length"].sum()  
    cov_dict[tr] = coverage

    temp_df = temp_df.sort_values(by="sstart").reset_index(drop=True)
    q_overlaps, q_cov = calc_overlaps(temp_df, "s") 
    right = q_overlaps["Right"].tolist()
    left = q_overlaps["Left"].tolist()
    distances = []
    if len(right)>1:
        for i in range(len(right)-1):
            distances.append(left[i+1]-right[i])
    else:
        distances = [0]
    dist_dict[tr] = distances

# all_blast.loc[all_blast["qseqid"] == tr, "Tr_coverage"] = coverage
all_blast["Tr_coverage"] = all_blast["qseqid"].apply(lambda x: cov_dict[x])
all_blast["s_distances"] = all_blast["qseqid"].apply(lambda x: dist_dict[x])
all_blast["max_distance"] = all_blast["s_distances"].apply(lambda x: max(x))
all_blast["srange"] = all_blast.apply(lambda x: [x.sstart, x.send], axis = 1)
all_blast["qrange"] = all_blast.apply(lambda x: [x.qstart, x.qend], axis = 1)
all_blast["counts"] = all_blast.groupby("qseqid")["sseqid"].transform("count")
```

# Step 7 - Sort dataframe by subject start

Self explanatory, performed for the sake of further processing. Run the cell.

```
cov50=all_blast.sort_values("sstart").reset_index(drop=True)
```

# Step 8 - Filter by maximal intron length and minimal transcript coverage

In this example, if the space between two neighbouring alignments exceeds 20000 bp, the transcript will be filtered out. Run the cell.

```
MAX_INTRON_LENGTH = 20000
cov50 = cov50[cov50["max_distance"]<MAX_INTRON_LENGTH]
```

![Screenshot](./images/large_intron_filter.png)

Furthermore, if the transcript is insufficiently covered by alignments (in this example, if less than 25%), it is filtered out. Program accounts for alignments that overlap. Run the cell.

```
MIN_COVERAGE = 25
cov50 = cov50[cov50["Tr_coverage"]>MIN_COVERAGE]
```

![Screenshot](./images/coverage_filter.png)

This step significantly reduces the number of false positives. Feel free to tweak MIN_INTRON_LENGTH and MIN_COVERAGE variables, based on your knowledge of organism biology.

Lastly, save the list of the transcript IDs remaining after filtration for further processing. Specify the save path in the FILTERED_TRANCRIPT_IDS variable.

```
FILTERED_TRANCRIPT_IDS = "./aligned_lrRNAseq_cov25_filtered.txt"
write_list(cov50.drop_duplicates("qseqid")["qseqid"].tolist(), FILTERED_TRANCRIPT_IDS)
```

# Step 9 - Extract the transcripts that passed the filtration and import into the notebook

Step produces a FASTA file including only the transcripts that passed the filtration procedure. Seqtk subseq extracts transcripts by transcript ID from the main lrRNAseq fasta. Enter path to the initial lrRNAseq file into the ALL_TRANSCRIPTS_PATH variable and the name of the fasta file that will contain post-filtered transcripts into the FILTERED_TRANSCRIPTS_FASTA_PATH variable. Run the cell.

```
ALL_TRANSCRIPTS_PATH = "./all_transcripts.fasta"
FILTERED_TRANSCRIPTS_FASTA_PATH = "filtered_transcripts.fasta"
os.system(f"seqtk subseq {ALL_TRANSCRIPTS_PATH} {FILTERED_TRANSCRIPT_IDS} > {FILTERED_TRANSCRIPTS_PATH}")
```

Import these transcripts into the notebook environment. If your transcripts still contain barcode sequences at either end (the example data contains barcodes of length 39 bp at either end), you need to remove the barcodes by specifying their length into the BARCODE_LEN variable. The program then imports these transcripts and calclates various statistics revolving around their Open Reading Frames.

```
BACRODE_LEN = 39
cov50_fasta = parse_fasta(FILTERED_TRANSCRIPT_FASTA_PATH,rm_barcode=BARCODE_LEN)[0]
all_stats = find_tr_ORF_stats(cov50_fasta["seq"].tolist(),cov50_fasta["id"].tolist(),tissue=True, include_seq=True)
```

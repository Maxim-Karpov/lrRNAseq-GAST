# Planned features

- In-depth view of splice junctions and general magnification - you will be able to select a transcript of choice from the main figure generated, and inspect closer, where the transcript's exons aligned to the genome.
- Visualisation of untranslated regions - similarly to exon visualisation, a different plot will be created displaying the location of 5' and 3' UTR alignments.
- Raw lrRNAseq processing protocol - a Jupyter-Notebook with a custom lrRNAseq barcode deconvolution protocol with superior performance than mainstream methods (e.g. PacBio Lima). Further processing of lrRNAseq to remove sequences of poor quality (e.g. no polyA) and orient the sequences in forward direction.
- Improved coverage filter - coverage filter which is based on just the coding regions.
- Improved ease of extracting transcripts of interest from the plot - ability to extract transcripts based on genomic coordinate range.
- Provide colour coding to the alignments based on alignment e-value, PID, or bitscore.
- Filter by e-value.

# Known bugs/design flaws

- Labelling of largest ORF/coding sequence may not be most accurate in the flanking exons.
- Optimal filtering protocol has not been determined yet. In case of finding difficult, highly variable, short exons, current filtering requires experimentation and understanding of biological sequences at hand as well as the software.

# Implemented

- lrRNAseq-GAST 1.0.

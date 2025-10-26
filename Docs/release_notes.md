# Planned features

- Alignment scoring systems based on sequence heuristics: Alignment scoring system based on similarity to splice-site consensus sequence. Short sequences flanking the alignment will be extracted and compared to a known splice-site consensus sequence for the organism to obtain a similarity score which will increase/decrease the existing e-value of the alignment. Alignment scoring based on presence/absence of start/stop codon at the beginning/end of alignment. Scoring based on intrinsic properties of the sequence e.g. kmer composition, GC content. Scoring system based on the structural sequentiality of the transcript i.e. alignments that sequentially start in 3'UTR move to coding region then to 5' UTR; coding regions scoring higher than UTRs. Heuristics based on the initial screening of the genomic sequence composition. Alignment depth as an indicator of confidence. Possible introduction of external data, such as protein sequences for model training/improvement. All of these would greatly improve the confidence of legitimate alignments, and allow for better filtering.
- Improved intron filter - an intron filter which will not remove the entire transcript alignment from the track, but rather split the transcript alignment at points of exceeding intron length into two (or more) separate alignments. This could potentially empower filtering process, and would help visually recognise legitimate transcript alignments.
- Automated re-blasting upon detection of a plausible gene region - e-value in BLAST is influenced by the length of the subject (genomic) sequence. When the genomic sequence is large, it could cause short and variable exons to go undetected as their e-value decreases. If upon detecting a plausible gene coding region, it is cut and re-BLASTed, the e-value will increase to make the exon become visible.
- In-depth view of splice junctions and general magnification - you will be able to select a transcript of choice from the main figure generated, and inspect closer, where the transcript's exons aligned to the genome.
- Visualisation of untranslated regions - similarly to exon visualisation, a different plot will be created displaying the location of 5' and 3' UTR alignments.
- Raw lrRNAseq processing protocol - a Jupyter-Notebook with a custom lrRNAseq barcode deconvolution protocol with superior performance than mainstream methods (e.g. PacBio Lima). Further processing of lrRNAseq to remove sequences of poor quality (e.g. no polyA) and orient the sequences in forward direction.
- Improved coverage filter - coverage filter which is based on just the coding regions.
- Improved ease of extracting transcripts of interest from the plot - ability to extract transcripts based on genomic coordinate range.
- Provide colour coding to the alignments based on alignment e-value, PID, or bitscore.
- Filter by e-value.
- User-friendly interface/wrapper.

# Known bugs/design flaws

- Labelling of largest ORF/coding sequence may not be most accurate in the flanking exons.
- Optimal filtering protocol has not been determined yet. In case of finding difficult, highly variable, short exons, current filtering requires experimentation and understanding of biological sequences at hand as well as the software.
- There is currently nothing implemented to indicate that the alignment happened in reverse orientation. Reverse alignments, in rare cases, can skew perception. Furthermore, the alignment region numbering system treats all types of alignments the same.

# Implemented

- lrRNAseq-TAST 1.0.

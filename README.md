# Response to Bean et al. (2018): Gain‐of‐function mutations in beet DODA2 identify key residues for betalain pigment evolution

Code in this repository is provided as-is and comes with absolutely no warranty. Feel free to share and use (at your own risk) but cite the paper when you do.

## Contents

- phylogenetics
  - data
    - **sequences.cds.fa**: coding sequences included in analysis
    - **sequences.pep.fa**: amino acid sequences included in analysis
    - **alignment.pep.fa.clustalo.aln**: clustalo alignment of amino acid sequences
    - **alignment.pep.fa.clustalo.aln-cln**: clustalo alignment of amino acid sequences, only >= 95% occupancy columns
    - **alignment.pep.fa.clustalo.aln-cln.treefile**: inferred tree, with 200 non-parametric bootstrap replicates
  - scripts
    - **calc_site_specific_divergence_aa.py**: script to calculate site-specific amino acid divergences
    - **parse_fasta.py**: utility
- chimeras
  - data
    - **chimera_data.csv**: CSV of chimera and fluorescence, expressed as fold-change over background (see analysis script for name correspondence)
  - scripts
    - **chimera_analysis.R**: R script for statistical analysis of chimera fluorescence data

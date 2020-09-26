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
- statistics
  - data
    - **fig_4a.csv**: CSV to reproduce fig 4a
    - **fig_4b.csv**: CSV to reproduce fig 4b
    - **fig_4c.csv**: CSV to reproduce fig 4c
    - **fig_5.csv**: CSV to reproduce fig 5. Note that name correspondences are given in statistics/scripts/analysis.R
  - scripts
    - **analysis.R**: R script for statistical analysis

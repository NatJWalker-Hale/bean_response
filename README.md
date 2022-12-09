# Response to Bean et al. (2018): Gain‐of‐function mutations in beet DODA2 identify key residues for betalain pigment evolution

Code in this repository is provided as-is and comes with absolutely no warranty. Feel free to share and use (at your own risk) but cite the paper when you do.

## Contents

- phylogenetics
  - data
    - **alignment.pep.fa.clustalo.aln**: clustalo alignment of amino acid sequences
    - **alignment.pep.fa.clustalo.aln-cln**: clustalo alignment of amino acid sequences, only >= 95% occupancy columns
    - **alignment.pep.fa.clustalo.aln-cln.treefile**: inferred tree, with 200 non-parametric bootstrap replicates
    - **alignment_BvDODAa2_corres.tsv**: correspondence between cleaned alignment positions and sites in BvDODAa2
    - **DODA1_DODA2_jsd.csv**: Jensen-Shannon Divergences between DODA1- and DODA2-like sequences in the cleaned alignment
    - **sequences.cds.fa**: coding sequences included in analysis
    - **sequences.pep.fa**: amino acid sequences included in analysis
  - scripts
    - **calc_site_specific_divergence_aa.py**: script to calculate site-specific amino acid divergences
    - **color_resi_pymol.py**: a PyMOL function to colour multiple residues based on values in a table
    - **DODA1_DODA2_colours_for_pymol.R**: script to get colours for PyMOL
    - **parse_fasta.py**: utility
- statistics
  - data
    - **Nbenthamiana_MN.csv**: CSV to reproduce Fig 4a
    - **Scerevisiae_GI.csv**: CSV to reproduce Fig 4b
    - **Scerevisiae_HC.csv**: CSV to reproduce Fig 4c
    - **DODA_activity_diff_pH.csv**: CSV to reproduce Fig 5a
    - **DODAa2_a2-mut3_activity_pH8.5.csv**: CSV to reproduce Fig 5b
    - **DODA_activity_pH6.csv**: CSV to reproduce Fig 5c
  - scripts
    - **analysis.R**: R script for statistical analysis
    - **curve_fitting.R**: R script for enzymatic analysis

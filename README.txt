# PeptideViz

PeptideViz provides the capability to visualize a prediction of how mutations will affect a set of
known peptides' binding affinity.

## System Requirements

python3
matplotlib
pandas
numpy

## Initial Data Requirements

The NGS sets of peptide data:
all_pep_seqs/All_peptides_Set1.csv
all_pep_seqs/All_peptides_Set2.csv
all_pep_seqs/All_peptides_Set3.csv

***If there are less than 3 sets, copy one of them to make Set3 (and Set2 if needed)

The set of known peptides:
python/og_peptides.txt

## Instructions

***Run all files from a command prompt

1. Run python/main_alg.py
    - This process will take hours, so it is ideal to let it run overnight

2. Run python/visualizer.py, python/visualizer_pos.py, and/or python/visualizer_3D.py
    - visualizer.py will display how a mutation from one amino acid to another affects binding
      affinity, regardless of the position of the mutation

    - visualizer_pos.py will display how a mutation to a given amino acid in a given position
      affects binding affinity, regardless of the initial amino acid in that position

    - visualizer_3D.py will display how a mutation from one amino acid to another in a given
      position affects binding affinity.

3. To visualize the standard deviation of the data from any of the visualizer files, navigate into
   the desired visualizer file and follow the instructions at the top of the file.

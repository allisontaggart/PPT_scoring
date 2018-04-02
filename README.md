# PPT_scoring
Polypyrimidine tract scoring script.

This script accepts fasta file of introns as input and identifies/scores polypyrimidine tracts

Script adapted from rules in: PMID11854178
Searches for PPT in 3' terminal 40 nucleotides
Nucleotides assigned values of -2 for A/G, +2 for C, +3 for U
Positional weighting of 1 for +/- 2 positions, 2 for +/- positions and 3 for position 0
Runs with >=2 are considered potential PPT
Edges pruned/extended by one nt to ensure that it starts/ends with pyrimidine
Runs of >9 nt are candidate PPT
Runs of 5-9 nt with >4 uridines are candidate PPT

Script outputs all candidate PPT in window and PPT scores

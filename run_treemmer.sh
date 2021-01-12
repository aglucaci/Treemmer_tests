#!/bin/bash

#git clone https://git.scicore.unibas.ch/TBRU/Treemmer.git
#mkdir -p testing
#cd testing

#1
#python ../Treemmer_v0.3.py testing/MACSE2_Out_Codons_recombinants_0.fas.treefile

#2
python ../Treemmer_v0.3.py testing/MACSE2_Out_Codons_recombinants_0.fas.treefile --verbose 1 -RTL 0.8

# end of file

#!/bin/bash

## @Usage parallel ::: "bash run_pipeline.sh"

for fasta in subsampled_fasta/*.fasta; do

    echo "# Processing: "$fasta

    ##does output exist?
     
    ## create a new alignment
    #java -jar macse_v2.05.jar -prog alignSequences -seq $fasta -gap_op -7 -gap_ext -1 -fs -30 -gc_def 1 -stop -100 -out_AA $fasta"_Out_AA.fas" -out_NT $fasta"_Out_Codons.fas"
    java -jar macse_v2.05.jar -prog alignSequences -seq $fasta -out_AA $fasta"_Out_AA.fas" -out_NT $fasta"_Out_Codons.fas"
    #java -jar macse_v2.05.jar -prog alignSequences -seq $fasta -out_AA $fasta"_Out_AA.fas" -out_NT $fasta"_Out_Codons.fas" -max_refine_iter 3 -local_realign_init 0.3 -local_realign_dec 0.2

    ## make a new tree
    FastTree -gtr -nt $fasta"_Out_Codons.fas" > $fasta"_Out_Codons.fas.fasttree"

    ##selection analyses
    #hyphy LIBPATH=$res fel --alignment $fasta"_Out_Codons.fas" --tree $fasta"_Out_Codons.fas.fasttree" --output $fasta"_Out_Codons.fas.FEL.json"
    hyphy fel --alignment $fasta"_Out_Codons.fas" --tree $fasta"_Out_Codons.fas.fasttree" --output $fasta"_Out_Codons.fas.FEL.json"
done 


# organize outputs.

mkdir -p macse2_realigned
mv *.fas macse2_realigned/

mkdir -p FastTree_nwks
mv *.fasttree FastTree_nwks/

mkdir -p FEL
mv *.FEL.json FEL/


#end of file

#!/bin/bash

#Declares
REFERENCE="p53/Human_p53.fasta"


#Main subroutine
for fasta in p53/subsampled_fasta/*.fasta; do
    echo "Processing: "$fasta

    hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $fasta --reference $REFERENCE --keep-reference No --remove-stop-codons Yes
    
    NUCFILE=$fasta"_nuc.fas"
    PROTFILE=$fasta"_protein.fas"
    
    mafft --auto $PROTFILE > $fasta"_protein.msa"
    
    CODONOUTPUT=$fasta"_codonmsa.fasta"
    hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $fasta"_protein.msa" --nucleotide-sequences $NUCFILE --output $CODONOUTPUT
    
    FastTree -gtr -nt $$CODONOUTPUT > $fasta"_codonmsa.nwk"
    
    hyphy fel --alignment $CODONOUTPUT --tree $fasta"_codonmsa.nwk" --output $CODONOUTPUT".FEL.json"

done

#Organize outputs
#mkdir -p NEWICKS
#mv p53/subsampled_fasta/*.nwk NEWICKS

mkdir -p p53/FEL
mv p53/subsampled_fasta/*.FEL.json p53/FEL


# End of file

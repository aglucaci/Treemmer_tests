#!/bin/bash

#Declares
REFERENCE="p53/Human_p53.fasta"


#Main subroutine
for fasta in p53/subsampled_fasta/*.fasta; do
    echo "# Processing: "$fasta

    NUCFILE=$fasta"_nuc.fas"
    PROTFILE=$fasta"_protein.fas"

    # Pre-MSA
    if [ -s $NUCFILE ]; then
        echo "# Nucleotide file exists"
    else
        hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $fasta --reference $REFERENCE --keep-reference No --remove-stop-codons Yes
    fi
    
    prunedNUCFILE=$fasta"_nuc_pruned.fas"
    prunedPROTFILE=$fasta"_protein_pruned.fas"
    
    # Pruned, checking for errors
    if [ -s $prunedNUCFILE ]; then
        echo "# Pruned file exists, checked for - and X's"
    else
        python Prune_nuc_and_protein.py $NUCFILE $PROTFILE
    fi
    
    # MAFFT
    if [ -s $fasta"_protein.msa" ]; then
        echo "# Protein MSA exists"
    else
        #mafft --auto $PROTFILE > $fasta"_protein.msa"
        mafft --auto $prunedPROTFILE > $fasta"_protein.msa"
    fi
    
    # Codon aware alignment
    CODONOUTPUT=$fasta"_codon.alignment"
    
    if [ -s $CODONOUTPUT ]; then
        echo "Codon alignment exists"
    else
        #hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $fasta"_protein.msa" --nucleotide-sequences $NUCFILE --output $CODONOUTPUT
        
        echo ""
        echo hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $fasta"_protein.msa" --nucleotide-sequences $prunedNUCFILE --output $CODONOUTPUT
        echo ""
        hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa $fasta"_protein.msa" --nucleotide-sequences $prunedNUCFILE --output $CODONOUTPUT
    fi
    
    # FastTree
    treeOUTPUT=$fasta".nwk"
    if [ -s $treeOUTPUT ]; then
        echo "FastTree exists"
    else
        echo FastTree -gtr -nt $CODONOUTPUT > $treeOUTPUT
        FastTree -gtr -nt $CODONOUTPUT > $treeOUTPUT
    fi
    
    felOUTPUT=$CODONOUTPUT".FEL.json"
    
    if [ -s $felOUTPUT ]; then
        echo "FEL output exists"
    else
        # Selection Analyses
        hyphy fel --alignment $CODONOUTPUT --tree $treeOUTPUT --output $CODONOUTPUT".FEL.json"
    fi

done

#Organize outputs
#mkdir -p NEWICKS
#mv p53/subsampled_fasta/*.nwk NEWICKS

mkdir -p p53/FEL
mv p53/subsampled_fasta/*.FEL.json p53/FEL


# End of file

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:04:15 2020

@author: alexander g. lucaci


The idea for this is that:
    
    Given a protein sequence and a transcript sequence
    I find the codons by stepping over the transcript sequence until the translated sequence matches the protein sequence
    that way, I have only the codons and not the additional sequences from the transcript
    (Which may be useful later)
    I will also create two output files
        One with the STOP codon stripped (this makes it hyphy compatible.)
        One with the STOP codons (may be useful later, codon bias?)
"""

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_rna
import os
import sys
from scipy import stats

# =============================================================================
# Declares
# =============================================================================
#PROTEIN = "human_protein_sequences.fasta"
#TRANSCRIPTS = "human_kinase_cDNA.fasta"
#OUTPUT = "human_kinase_cDNA_codons.fasta"

PROTEIN = sys.argv[1]
TRANSCRIPTS = sys.argv[2]
OUTPUT = sys.argv[3]
ERRORS = {}

# =============================================================================
# Helper functions
# =============================================================================
#turn into class

def Process(PROTEIN_DESC, PROTEIN_SEQ, TRANSCRIPT_DESC, TRANSCRIPT_SEQ):
    #print([PROTEIN_DESC], [TRANSCRIPT_DESC])
    
    # print()
    #print(PROTEIN_DESC, "\nProtein seq length (AA):", len(PROTEIN_SEQ))
    # print(PROTEIN_DESC, PROTEIN_SEQ)
    #print()
    
    #Loop over TRANSCRIPT_SEQ
    start = 0
    
    NT_SEQ_LENGTH = len(PROTEIN_SEQ) * 3
    
    while start < len(str(TRANSCRIPT_SEQ)):
        try:
            #[including: up to but excluding]
            coding_dna = TRANSCRIPT_SEQ[start : start + NT_SEQ_LENGTH] #translated
            if len(str(coding_dna)) % 3 == 0:
                coding_dna = TRANSCRIPT_SEQ[start : start + NT_SEQ_LENGTH].translate(table='Standard')
            else:
                start += 1
                continue
            #print(coding_dna)
        except:
            pass
        
        #if start == 202:
        #    print(start, coding_dna, "\n", len(coding_dna))
        #print(coding_dna)
        if coding_dna == str(PROTEIN_SEQ):
            #print("\n#### FOUND", coding_dna)
            #print("#### CODONS", TRANSCRIPT_SEQ[start: start + NT_SEQ_LENGTH + 3]) # has stop codon
            #print("\n#### CODONS", TRANSCRIPT_SEQ[start: start + NT_SEQ_LENGTH]) # NO stop codon
            break
        else:
            pass
            # print()
        
        start += 1
        #if start == 301: break
        #end if
    #end while
    
    return TRANSCRIPT_SEQ[start: start + NT_SEQ_LENGTH]
#end method

# =============================================================================
# Main subroutine.
# =============================================================================
def progressBar(value, endvalue, bar_length=75):
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()
#end method


def main(): # Really to verify things.
    global PROTEIN, TRANSCRIPTS
    print("\tTRANSCRIPT INPUT FILE:", TRANSCRIPTS)
    print("\tPROTEIN INPUT FILE:", PROTEIN)
    print()
    
    protein_list = []
    transcript_list = []
    
    with open(TRANSCRIPTS, "r") as handle:
        #x = SeqIO.parse(handle, "fasta")
        #print(len(x))
        trans_count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            trans_count +=1
            transcript_list.append(record.description)
        print("\tTranscripts:", trans_count)    
    handle.close()
    
    with open(PROTEIN, "r") as handle:
        #x = SeqIO.parse(handle, "fasta")
        #print(len(x))
        prot_count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            prot_count +=1
            protein_list.append(record.description)
        print("\tProteins:", prot_count)
    handle.close()
    
    #assert(trans_count == prot_count, "Counts do not match") # Check to make sure we have the same number of proteins and transcripts
    return trans_count, prot_count
    #for n, item in enumerate(transcript_list):
    #    print(item, [protein_list[n]])
#end method

# =============================================================================
# Main
# =============================================================================
#Verify files exist
print("# =============================================================================")
print("# Processing... ")
trans_count, prot_count = main()
print("# =============================================================================")
#Looks like species all match up in transcript and protein fasta.
#This is exceptional, will need to look for species name (from protein desc.) in transcript desc.

# DEBUG
# sys.exit(1)
      
# Create empty output file.
print("# Saving output to:", OUTPUT)
with open(OUTPUT, "w") as fh:
    fh.write("")
fh.close()


# Main program
successful_count = 0 
num_errors = 0
errors_IDs = []

print("# Opening", PROTEIN, "file")
with open(PROTEIN, "r") as prot_handle:
    for n, record in enumerate(SeqIO.parse(prot_handle, "fasta")):
        
        #if n == 1: break
        
        #print("\n" + str(n))
        
        # Grab protein data.
        protein_id = record.id 
        protein_desc = record.description
        protein_seq = record.seq
        
        # print("# Opening", TRANSCRIPTS, "file")
        with open(TRANSCRIPTS, "r") as transcript_handle:
            for m, transcript_record in enumerate(SeqIO.parse(transcript_handle, "fasta")):
                if m == n:
                    progressBar(n, trans_count)
                    
                    # Grab Transcript Data
                    transcript_id = transcript_record.id
                    transcript_desc = transcript_record.description
                    transcript_seq = transcript_record.seq
                #end if
            #end inner for
        transcript_handle.close()
        #end inner with
        
        # Process
        #assert(str(protein_desc) == str(transcript_desc), "Does not match")
        #print(protein_desc, [transcript_desc])
        codons = ""
        codons = Process(protein_desc, protein_seq, transcript_desc, transcript_seq)
        
        #if "XM_017345" in transcript_desc: print("\n", transcript_desc, [codons], len(codons), str(codons))
        #assert(len(codons) > 0), "EMPTY CODONS: " + protein_desc
        if len(codons) == 0: #ERROR
            num_errors += 1
            errors_IDs += [transcript_desc]
            ERRORS[transcript_desc] = {"Protein": protein_seq, "Transcript": transcript_seq}
            continue
        
        successful_count += 1
        #Print out transcript desc and TRIMMED codons transcript.
        with open(OUTPUT, "a") as fh:
            fh.write(">" + transcript_desc + "\n" + str(codons) + "\n")
        fh.close()
        
    #end outer for
prot_handle.close()


#end outer with

print("\n# Found:", successful_count)

print()
print("## Errors ##", num_errors)

for item in errors_IDs: print(item)

print()

print("## Detailed Errors printout ##")
for item in ERRORS.keys():
    print(item)
    print("Protein:", ERRORS[item]["Protein"])
    if "X" in ERRORS[item]["Protein"]: 
        #print [pos for pos, char in enumerate(s) if char == c]
        X_chars = [pos for pos, char in enumerate(ERRORS[item]["Protein"]) if char == "X"]
        print("## ERROR #1: There is an 'X' character in the protein sequence.", "At position(s):", X_chars)
        
        
    print("Transcript:", ERRORS[item]["Transcript"])
    print()

# =============================================================================
# End of file    
# =============================================================================

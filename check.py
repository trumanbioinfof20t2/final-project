'''
A program that checks for validity of sequences, with ambiguity codes, and 
then aligns sequences and calculates a distance matrix
Last Modified: 11/17/20
'''

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import sys

from readfasta import readfasta

'''
checksequence - a function that checks if a sequence only contains valid amino acids
@param sequence_input - the sequence which we are performing the check on
@return True, if the AA sequence is valid, false otherwise
'''
def checksequence(sequence_input):
    # valid amino acids are all these letters http://130.88.97.239/bioactivity/aacodefrm.html
    validAA = 'GALMFWKQESPVICYHRNDT-'
    
    #sequence = input('Please enter a sequence: ').upper()
    sequence = sequence_input.upper()
    
    # will print True if every character in the sequence belongs to valid_dna
    if not all(i in validAA for i in sequence):
        return True
    else:
        return False
'''
align - a function that takes an unaligned FASTA file and runs it through clustal, 
    creating a aligned clustal file
@param unalignedFileName - the file name of the unaligned FASTA file of sequences
'''
def clustalAlign(unalignedFileName):
    cmd = ClustalwCommandline("clustalw2", infile=unalignedFileName)
    print(cmd)

    # this line sends the clustal command to the terminal
    stdout, stderr = cmd()

'''
distanceCalculation - a function that uses Biopython's in-built distance calculator, and saves
    the distance matrix to a txt file
@param alignment - a clustal file that has all our sequences aligned
@return the distance matrix that will be used in tree construction
'''
def distanceCalculation(alignment):
    distanceFileName = "Distance.txt"

    # can also use 'blosum62' instead of identity
    calculator = DistanceCalculator('identity')
    distanceMatrix = calculator.get_distance(alignment)

    sys.stdout = open(distanceFileName, "w")
    print(distanceMatrix)
    sys.stdout.close()

    return distanceMatrix

def main():
    
    #just some filler variables for file names
    unalignedFileName = "sequences_aa.fa"
    alignedFileName = "sequences_aa.aln"

    #can ask user for input here at final version
    sequences = readfasta(unalignedFileName)

    #checks the seqeunces for validty with amino acid letters
    max_len = 0;
    for sequence in sequences:
        if len(sequence[1]) > max_len:
            max_len = len(sequence[1])
        if checksequence(sequence[1]) == True:
            print(sequence[0], "is not valid")
        else:
            print(sequence[0], "is valid!")

    clustalAlign(unalignedFileName)

    # prints the alignment file to the terminal
    alignment = AlignIO.read(alignedFileName, "clustal")
    print(alignment)

    distanceMatrix = distanceCalculation(alignment)

if __name__ == "__main__":
    main()
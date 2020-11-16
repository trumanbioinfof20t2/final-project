'''
A program that checks for validity of sequences, with ambiguity codes, and 
then aligns sequences and calculates a distance matrix

'''

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import sys


from readfasta import readfasta

def checksequence(sequence_input):
    # valid DNA nucleotides A, C, G, and T.
    valid_dna = 'ACGTRYSWKMBDHVN-'
    #sequence = input('Please enter a sequence: ').upper()
    sequence = sequence_input.upper()
    # will print True if every character in the sequence belongs to valid_dna
    if not all(i in valid_dna for i in sequence):
        return True
    else:
        return False

def main():
    '''
    Reads a FASTA file,
    '''
    # Ask for a file name and read in the FASTA file of sequences
    #fileName = input("Please enter a FASTA file name: ")
    fileName = "sequences.fa"
    sequences = readfasta(fileName)

    max_len = 0;
    for sequence in sequences:
        if len(sequence[1]) > max_len:
            max_len = len(sequence[1])
        if checksequence(sequence[1]) == True:
            print(sequence[0], "is not valid")
        else:
            print(sequence[0], "is valid!")
    print("sequences are valid!")

    cmd = ClustalwCommandline("clustalw2", infile=fileName)
    print(cmd)

    #stdout, stderr = cmd()

    align = AlignIO.read("sequences.aln", "clustal")
    print(align)

     # can also use 'blosum62' instead of identity
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(align)

    sys.stdout = open("Distance.txt", "w")
    print(dm)
    sys.stdout.close()

if __name__ == "__main__":
    main()

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Align

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
    sequences = readfasta("sequences.fa")

    max_len = 0;
    for sequence in sequences:
        if sequence[1].length > max_len
            max_len = sequence[1].length
        if checksequence(sequence[1]) == True:
            print(sequence[0], "is not valid")
        else
            print("rest of the sequences are valid!")
    print("sequences are valid!")

    for sequence in sequences:
        if sequence[1].length < max_len
            need_len = max_len - sequences[1].length
            sequence.append("-"*need_len)

    



        


    # aligner = Align.PairwiseAligner()
    # alignments = aligner.align("ATCTT", "ACTGAAT")

    # for alignment in sorted(alignments):
    #     print("Score = %.1f:" % alignment.score)
    #     print(alignment)



    # alignment = AlignIO.parse(open("sequences.fa"), "fasta")
    # print("Alignment length %i" % alignment.get_alignment_length())
    # for record in alignment:
    #     print(record.seq + " " + record.id)



    # aln = AlignIO.read(open(''), 'phylip')
    # print(aln)



if __name__ == "__main__":
    main()

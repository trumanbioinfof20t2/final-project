'''
Bioinformatics Team 2 Final Project
Miller, Ingli, Shahi, Winistoerfer
Last Modified: 11/17/2020
1. Check validity of amino acid sequences
2. Align sequences with clustal
3. Calculate a distance matrix
4. Create a series of phylogenetic trees
'''

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
from Bio.Phylo import draw as phyloDraw
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import matplotlib.pyplot as plt
import sys
import time

from readfasta import readfasta

'''
isValidSequence - a function that checks if a sequence only contains valid amino acids based on the single-character code
@param sequence_input - the sequence which we are performing the check on
@return True if the AA sequence is valid, False otherwise
'''
def isValidSequence(sequence_input):
    # valid amino acids are all these letters http://130.88.97.239/bioactivity/aacodefrm.html
    validAA = 'GALMFWKQESPVICYHRNDT-X'
    
    #sequence = input('Please enter a sequence: ').upper()
    sequence = sequence_input.upper()
    
    # return True if all characters in the sequence are valid
    # return False if at least one character is not valid
    return all(i in validAA for i in sequence)

'''
clustalAlign - a function that takes an unaligned FASTA file and runs it through clustal, 
    creating an aligned clustal file. Clustal writes this to a file with 
    the same name as the infile, but with the .aln extension.
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
@param distanceMethod - the matrix for the DistanceCalculator to use. Default: identity
@return the distance matrix that will be used in tree construction
'''
def distanceCalculation(alignment, distanceMethod='identity'):
    calculator = DistanceCalculator(distanceMethod)
    distanceMatrix = calculator.get_distance(alignment)

    return distanceMatrix

'''
buildNeighborJoiningTree - a function that uses Biopython's DistanceTreeConstructor to 
    build a representation of a Neighbor Joining tree based on the given distance matrix. 
    Sequence information is contained within the matrix object.
@param distanceMatrix - the calculated distance matrix for a set of aligned sequences.
@see distanceCalculation()
@return Tree object representing the NJ tree
'''
def buildNeighborJoiningTree(distanceMatrix):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distanceMatrix)
    tree.ladderize()
    return tree

'''
buildUPGMATree - a function that uses Biopython's DistanceTreeConstructor to 
    build a representation of a UPGMA tree based on the given distance matrix. 
    Sequence information is contained within the matrix object.
@param distanceMatrix - the calculated distance matrix for a set of aligned sequences.
@see distanceCalculation()
@return Tree object representing the UPGMA tree
'''
def buildUPGMATree(distanceMatrix):
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distanceMatrix)
    tree.ladderize()
    return tree

'''
buildParsimonyNNITree - a function that uses Biopython's ParsimonyTreeConstructor, 
    ParsimonyScorer, and NNITreeSearcher to build a representation of a Maximum 
    Parsimony tree based on the given alignment, starting tree, and optional 
    scoring matrix.
@param alignment - A sequence alignment in the appropriate Biopython object
@param startingTree - a Tree object to start NNI with, usually generated from UPGMA or NJ
@param scoreMatrix - Parsimony Scoring Matrix (in _Matrix object form) to use in the 
    Sankoff algorithm. Default: None (use the Fitch algorithm instead)
'''
def buildParsimonyNNITree(alignment, startingTree, scoreMatrix=None):
    scorer = ParsimonyScorer(scoreMatrix)
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher, startingTree)
    tree = constructor.build_tree(alignment)
    return tree

'''
drawTree - a function that draws a pictoral representation of a phylogenetic tree 
    using matplotlib and Biopython, then saves the image as a png file.
@param tree - Tree object to draw
@param title - Title of the tree to display
@param fileName - file name to save tree image to. Must be a .png
'''
def drawTree(tree, title, fileName):
    fig = plt.figure(figsize=(20,10), dpi=100)
    plt.suptitle(title, fontsize=40.0)
    axes = fig.add_subplot(1,1,1)
    plt.tight_layout()
    phyloDraw(tree, axes=axes, branch_labels=lambda c: round(c.branch_length, 3), do_show=False)
    plt.savefig(fileName)

# ------------------------------------------------------------------------------------------------------

def main():

    # Record starting time
    startTime = time.time()
    
    # Key File Names
    unalignedFileName = "sequences_aa2.fa"
    alignedFileName = "sequences_aa2.aln"
    upgmaTreeFileName = "upgma.png"
    njTreeFileName = "nj.png"
    parsimonyTreeFileName = "nniparsimony.png"

    # Read in the unaligned amino acid sequences FASTA file, single-letter code
    sequences = readfasta(unalignedFileName)

    # -------------------------------------

    # Step 1: Check the sequences for validity with amino acid single-letter code

    print("Checking sequences...")

    allSeqsValid = True
    for sequence in sequences:
        if not isValidSequence(sequence[1]):
            print(sequence[0], "is not valid.")
            allSeqsValid = False
    if not allSeqsValid:
        sys.exit("Invalid sequence(s) detected")

    print("All sequences valid.")

    # -------------------------------------

    # Step 2: Use Clustal to align the sequences

    print("Aligning Sequences using Clustal...")
    clustalAlign(unalignedFileName)

    # Read clustal alignment
    print("Loading Clustal alignment from disk...")
    alignment = AlignIO.read(alignedFileName, "clustal")

    print("Alignment Complete.")

    # -------------------------------------

    # Step 3: Calculate the distance matrix
    # first with the identity distance matrix
    print("Calculating identity distance matrix...")
    distanceMatrixIdentity = distanceCalculation(alignment, "identity")

    # then with the blosum62 distance matrix
    print("Calculating BLOSUM62 distance matrix...")
    distanceMatrixBlosum62 = distanceCalculation(alignment, "blosum62")

    print("Distance matrices calculated.")

    # -------------------------------------

    # Step 4: Create a series of phylogenetic trees

    # Step 4.1: Build UPGMA Trees
    # First with the identity matrix
    print("Building UPGMA Identity tree...")
    upgmaIdentityTree = buildUPGMATree(distanceMatrixIdentity)
    drawTree(upgmaIdentityTree, "UPGMA - IDENTITY", "identity-"+upgmaTreeFileName)

    # Then with the BLOSUM62 matrix
    print("Building UPGMA BLOSUM62 tree...")
    upgmaBlosum62Tree = buildUPGMATree(distanceMatrixBlosum62)
    drawTree(upgmaBlosum62Tree, "UPGMA - BLOSUM62", "BLOSUM62-"+upgmaTreeFileName)

    # ---------------

    # Step 4.2: Build Neighbor Joining Trees
    # First with the identity matrix
    print("Building Neighbor Joining Identity tree...")
    njIdentityTree = buildNeighborJoiningTree(distanceMatrixIdentity)
    drawTree(njIdentityTree, "Neighbor Joining - IDENTITY", "identity-"+njTreeFileName)

    # Then with the BLOSUM62 matrix
    print("Building Neighbor Joining BLOSUM62 tree...")
    njBlosum62Tree = buildNeighborJoiningTree(distanceMatrixBlosum62)
    drawTree(njBlosum62Tree, "Neighbor Joining - BLOSUM62", "BLOSUM62-"+njTreeFileName)

    # ---------------

    # Step 4.3: Build Parsimony Trees using Nearest Neighbor Interchange
    
    # Step 4.3.1: With the Fitch Algorithm (no scoring matrix)

    # Step 4.3.1.1: Use UPGMA - Identity as a starting tree
    print("Building NNI Parsimony tree using Fitch Algorithm starting with UPGMA Identity...")
    parsimonyUPGMAIdent = buildParsimonyNNITree(alignment, upgmaIdentityTree)
    drawTree(parsimonyUPGMAIdent, "Parsimony NNI - Start with UPGMA IDENTITY", "upgmaident-"+parsimonyTreeFileName)

    # Step 4.3.1.2: Use UPGMA - BLOSUM62 as a starting tree
    print("Building NNI Parsimony tree using Fitch Algorithm starting with UPGMA BLOSUM62...")
    parsimonyUPGMABlosum62 = buildParsimonyNNITree(alignment, upgmaBlosum62Tree)
    drawTree(parsimonyUPGMABlosum62, "Parsimony NNI - Start with UPGMA BLOSUM62", "upgmablosum62-"+parsimonyTreeFileName)

    # Step 4.3.1.3: Use Neighbor Joining - Identity as a starting tree
    print("Building NNI Parsimony tree using Fitch Algorithm starting with NJ Identity...")
    parsimonyNJIdent = buildParsimonyNNITree(alignment, njIdentityTree)
    drawTree(parsimonyNJIdent, "Parsimony NNI - Start with NJ IDENTITY", "njident-"+parsimonyTreeFileName)

    # Step 4.3.1.4: Use Neighbor Joining - BLOSUM62 as a starting tree
    print("Building NNI Parsimony tree using Fitch Algorithm starting with NJ BLOSUM...")
    parsimonyNJBlosum62 = buildParsimonyNNITree(alignment, njBlosum62Tree)
    drawTree(parsimonyNJBlosum62, "Parsimony NNI - Start with NJ BLOSUM62", "njblosum62-"+parsimonyTreeFileName)

    # -------

    # Step 4.3.2: With the Sankoff Algorithm and the X matrix
    # TODO

    print("Tree building complete.")
    # -------------------------------------

    # Record end time
    endTime = time.time()
    elapsedSeconds = endTime - startTime
    print("Completed in {:.3f} seconds".format(elapsedSeconds))

# If this file is executed, run main()
# Don't run main() if this file is imported elsewhere
if __name__ == "__main__":
    main()
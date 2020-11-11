'''
Jon Beck
A function to read a fasta file
Last modified: 26 August 2020
'''

'''
readfasta - a function to read a fasta file
Parameter: a filename that must be in fasta format.  The file is
assumed to have:
1. the first line a header line
2. arbitrary blank lines
3. every line (especially including the last) is terminated by a newline
   terminator (enter key)
4. no line has only spaces on it

Return: a list of lists. Each inner list will have two elements:
1. the sequence identifier, i.e., the characters between the leading ">"
   and the first space
2. the sequence, a single string of all the letters with no line terminators
'''
def readfasta(filename):
    result_list = []
    with open(filename, 'r') as infile:

        # process the first line, which must be a header line
        line = infile.readline()
        header_line = line.rstrip()
        label = header_line[1:].split(' ')[0]

        # initialize the sequence accumulator
        sequence = ''

        # process all the rest of the lines in the file
        for line in infile:
            line = line.rstrip()

            # ignore blank lines
            if line != '':

                # if it's a header line, finish the previous sequence
                # and start a new one
                if line[0] == '>':
                    result_list.append([label, sequence])

                    label = line[1:].split(' ')[0]
                    sequence = ''
            
                    # if we're here, we must be in letters of the sequence
                else:
                    sequence += line
            
    # we're done, so clean up, terminate the last sequence, and return
    result_list.append([label, sequence])
    return result_list


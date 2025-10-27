def fasta_load(filename):

"""
Loads fasta file as list of sequences. Also provides a list of header indices that could be returned as a variable should the user wish to do so.
"""
    try:
        # Open file
        with open(filename) as f:
            # Initialize list for sequences and read in file
            seqlist = []
            lines = f.readlines()
            # Initialize list for header indices.
            headerindices = []
            # Iterate through lines in fasta file, append stripped lines to seqlist and header lines to headerindices list.
            for line in lines:
                seqlist.append(line.strip())
                if line.startswith(">") == True:
                    headerindices.append(lines.index(line))
                else:
                    continue
            # Initialize final sequences list.
            sequences = []
            # Iterate through the headerindices. Grab all sequences between headers and join them. This only works as there is only one copy of the                       sequence between headers and nothing else.
            for i in range(len(headerindices)):
                if i < len(headerindices) - 1:
                    j = headerindices[i+1]
                else:
                    j = len(seqlist)
                x = "".join(seqlist[(headerindices[i]+1):j])
                sequences.append(x)
        print(sequences)
        return sequences
            
    except FileNotFoundError:
        print("Cannot find the files")

def validate_sequence(sequence):
    """ 
    Validates sequence by sum-checking: if all valid nucleotides, sum of nucleotides should equal the length of the sequence.
    """
    # Iterates over every sequence given. Checks for equality between nucleotide sum and length.
    for seq in sequence:
        if (seq.count("A") + seq.count("T") + seq.count("C") + seq.count("G")) == len(seq):
            print("The sequence contains only nucleotides A,T,G,C.")
        else:
            print("Provide a valid sequence")

def Tm_calculation(sequence):

"""
Calculates Tm values for given sequences.
"""
    # Initialize Tmlist
    Tmlist = []
    # Iterate over every sequence and count A, T, C, G occurences.
    for seq in sequence:
        w = seq.count("A")
        x = seq.count("T")
        y = seq.count("G")
        z = seq.count("C")
        # Simplified formula for Tm calculation.
        Tm = 64.9 + (41*((y+z-16.4)/(w+x+y+z)))
        Tmlist.append(Tm)
        print(Tmlist)
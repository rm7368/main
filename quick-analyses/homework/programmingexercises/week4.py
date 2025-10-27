def myfunc(file):
    with open(file) as f:
        fastadict = {}
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
                fastadict[header] = ""
            elif header:
                fastadict[header] += line 
    return(fastadict)

def validate_base_seq(nt_seq):
    seq = nt_seq.upper()
    test_logic = len(seq) == seq.count('A') + seq.count('C') + seq.count('T') + seq.count('G')
    try:
        assert test_logic, "not dna"
    except AssertionError:
        return("This is not a valid DNA sequence.")
    

fastadict = myfunc("programmingexercises/Bogus.fa")
validate_base_seq(fastadict["GeneA|Bogus Gene Name"])
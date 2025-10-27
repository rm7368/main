import argparse

# Function to parse the fastq file
def parse_fastq(filename):
    """
    Parses fastq sequence file, returns list of dictionaries. Each dictionary contains id, sequence, and quality entries.
    """
    with open(filename, "r") as f:
        # Open file, read lines and initialize list and dictionary.
        lines = f.readlines()
        dictlist = []
        seqdict = {}
        # Clean up lines, removing trailing new line char
        for i in lines:
            new = i.strip("\n")
            lines[lines.index(i)] = new
        # Find sequence id line and verify that it has a corresponding spacer line.
        for i in lines:
            j = lines.index(i) + 2
            k = lines.index(i)
            if i.startswith("@") and lines[j].startswith("+"):
    
                # Add id, sequence, and quality entries, then append to list            
                seqdict["id"] = lines[k]
                seqdict["sequence"] = lines[k+1]
                seqdict["quality"] = lines[j+1]
                dictlist.append(seqdict)
                seqdict = {}
            else:
                continue
        print(dictlist)
        return(dictlist)

def calculate_average_quality(quality_str):
    """
    Calculates the average quality score given a quality string from a FASTQ file, in an ASCII format.
    """
    # Initialize score list.
    charscorelist = []
    # Convert ASCII to Phred score and append to list.
    for char in quality_str:
        phred = ord(char) - 33
        charscorelist.append(phred)
    # Calculate average of Phred score list.
    avg = sum(charscorelist)/len(charscorelist)
    print(f"The average quality score of this string is {avg}")
    return avg

def filter_sequences(sequences, threshold):
    """
    Filter bad sequences out of a parsed list of fastq sequences.
    """
    # Initialize list of sequences.
    goodseqs = []
    # Calculate quality score, append to list if over threshold.
    for i in sequences:
        if calculate_average_quality(i["quality"]) >= threshold:
            goodseqs.append(i)
        else:
            continue
    print(goodseqs)
    return(goodseqs)

def write_fastq(sequences, output_filename):
    """
    Writes filtered fastq sequences to a new fastq file.
    """
    # Creates file, write each line individually for each sequence in sequences list.
    with open(output_filename, "w") as file:
        for i in sequences:
            file.write(i["id"] + "\n")
            file.write(i["sequence"] + "\n")
            file.write("+" + i["id"].strip("@") + "\n")
            file.write(i["quality"] + "\n")
            
def main():

    """
    Defines argument parser for function.
    """
    # Initialize argument parser and add arguments.
    parser = argparse.ArgumentParser(description="A script that filters a list of FASTQ sequences based on their Phred quality scores.")
    parser.add_argument("--input", type=str, help="Path to the input file.")
    parser.add_argument("--output", type=str, help="Path to the output file.")
    parser.add_argument("--threshold", type=int, help="Threshold score. All sequences with a higher average Phred quality score will be written to the new file.")
    args = parser.parse_args()

    dictlist = parse_fastq(args.input)
    filtered = filter_sequences(dictlist, args.threshold)
    write_fastq(filtered, args.output)
    
    print(f"Total sequences: {len(dictlist)}")
    print(f"{(((len("filtered_sequencing_results.fastq"))-1)/4)} sequences passed filtering.")
    print(f"{len(dictlist)-(((len("filtered_sequencing_results.fastq"))-1)/4)} sequences did not pass filtering.")

# Only run when conducted as command line.
if __name__ == "__main__":
    main()
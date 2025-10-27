# Initialize empty dictionary
locusdict = {}

# Open file with "with" statement
with open("/home/rm7368/programmingexercises/Athaliana_1000_annotation.txt") as file:

    # Read lines in, save header if needed
    lines = file.readlines()
    header = lines[0]

    # Read lines after header, split by tab, set variable values
    for i in range(1,1000):
        line = str(lines[i])
        line = line.split("\t")
        locus = line[1]
        besthit = line[10]

        # If the besthit value is empty, pass. If it is not, add to dict and continue. 
        if str(besthit) == "":
            pass
        else:
            locusdict[locus] = besthit
            continue

# Print locusdict
print(locusdict)

# Print locusdict accessed at specific locus.
print(locusdict["AT1G01220"])
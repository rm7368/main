locusdict = {}
with open("/home/rm7368/programmingexercises/Athaliana_1000_annotation.txt") as file:
    lines = file.readlines()
    header = lines[0]
    for i in range(1,1000):
        line = lines[i]
        line = str(line)
        line = line.strip().split("\t")
        print(line[10])
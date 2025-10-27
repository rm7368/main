from Bio import Entrez
import sys

Entrez.email = "Razvan.Matei@nyulangone.org"

try:
    print("Testing NCBI connection...")
    handle = Entrez.efetch(db="nucleotide", id="NM_010275", rettype="fasta", retmode="text")
    result = handle.read()
    handle.close()
    print("SUCCESS: NCBI connection works")
    print(f"Retrieved {len(result)} characters")
except Exception as e:
    print(f"FAILED: {e}")
    sys.exit(1)
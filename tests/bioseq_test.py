from Bio.Seq import Seq
from Bio.Restriction import *
from Bio.Restriction.Restriction import  RestrictionType

# Create a DNA sequence as a Seq object
dna_sequence = Seq("GAATTCAGTCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG")

class MyEnzyme(RestrictionType):
    """A custom restriction enzyme class."""
    # The name of the enzyme
    name = "MyEnzyme"
    # The recognition site of the enzyme
    site = "GAATTC"
    # The cut position of the enzyme
    cut_pos = 1

# Create a restriction batch with the custom enzyme
rb = RestrictionBatch([MyEnzyme()])


# Analyze the DNA sequence with the restriction batch
results = rb.search(dna_sequence)

# Print the results
for enzyme, cut_sites in results.items():
    print(f"{enzyme}: {cut_sites}")
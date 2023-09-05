"""
This script will algorithmically generate genomic edits to be applied to a reference genome in order to create edited
genomes for validation testing. The edits will be "recorded" in a CSV file before they are applied. The script will
be able to randomly generate edits within strict parameters specified by the user. Fields in the CSV that need to
generated include:

1) ID - Arbitrary identifying name

2) Edit type - Substitution, Deletion, Insertion, (consider ability to create SNPs for future use)

3) Start - Position in genome where edit sequence starts
    - TODO: investigate the fields in the inscripta CSV (how much flanking sequence is included?)

4) Stop -  Position in genome where edit sequence ends

5) Ref sequence - original sequence in reference genome between the start and the stop positions.

6) Edit sequence - sequence that replaces the reference sequence between the start and stop positions.

Besides these fields, the edit parameters need to be strictly defined:

a) Length of edit - range of length values the edits can take

b) Distribution - vague, but we need to decide how the edits are allowed to be distributed. No edits within
    n base pairs of the start or end of chromosome, for example.

c) Types of edits - Substitution, Deletion, or Insertion. Frequency for each can be a weight selected. For testing we
    can make all of them equally likely.

d) Number of edits - the user can specify how many edits should be generated.
"""


import csv
import random
from Bio import SeqIO


# Function to load the reference genome from a FASTA file
def load_reference_genome(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return record.seq, len(record.seq)


# Define edit types and their weights
edit_types = ["Substitution", "Deletion", "Insertion"]
edit_weights = [1, 1, 1]  # Equal weights for testing


def generate_substitution(start_pos, stop_pos, length):
    """
    This function generates an edit of the substitution type. The reference sequence in the edit range will be
    replaced by a different sequence that for testing purposes will be a sequence of "N"s. Will take
    edit parameters as arguments and use these to create the edit.
    """
    edit_sequence = "N" * length
    return edit_sequence


def generate_deletion(start_pos, stop_pos, length):
    """
    This function generates a deletion edit. The reference sequence in the edit range will be deleted. Will take
    edit parameters as arguments and use these to create the edit.
    """
    return ""


def generate_insertion(start_pos, stop_pos, length):
    """
    This function generates an insertion edit. After the first base at the start position, a new edit sequence will
    be inserted. For testing, the sequence will consist of "N"s. Will take
    edit parameters as arguments and use these to create the edit.
    """
    edit_sequence = "N" * length
    return edit_sequence


# Function to generate a random edit record
def generate_edit(chromosome_length):
    edit = {
        "ID": "ID-XXXX",  # You can generate a unique ID
        "Edit type": random.choices(edit_types, weights=edit_weights)[0],
        "Start": random.randint(100, chromosome_length - 100),
    }
    # Generate a random length within a specific range
    edit_length = random.randint(10, 50)  # Adjust the range as needed
    # Calculate the stop position based on the start and length
    edit["Stop"] = edit["Start"] + edit_length
    return edit


# Number of edits to generate (user-defined)
num_edits = 50

# Input FASTA file (provide the actual path, eventually this will be a user-defined argument)
fasta_file = "/home/casper/PycharmProjects/fantastic-lamp/validation testing/reference genome/lambda_phage.fasta"

# Load reference genome and its length
reference_sequence, chromosome_length = load_reference_genome(fasta_file)

# Generate and store edit records
edits = []
for _ in range(num_edits):
    edit = generate_edit(chromosome_length)
    start = edit["Start"]
    stop = edit["Stop"]

    if edit["Edit type"] == "Substitution":
        edit_sequence = generate_substitution(start, stop, stop - start + 1)
    elif edit["Edit type"] == "Deletion":
        edit_sequence = generate_deletion(start, stop, stop - start + 1)
    elif edit["Edit type"] == "Insertion":
        edit_sequence = generate_insertion(start, stop, stop - start + 1)

    edit["Ref sequence"] = str(reference_sequence[start - 1: stop])  # Adjust for 0-based indexing
    edit["Edit sequence"] = edit_sequence
    edits.append(edit)

# Write edits to a CSV file
with open("genomic_edits.csv", "w", newline="") as csvfile:
    fieldnames = ["ID", "Edit type", "Start", "Stop", "Ref sequence", "Edit sequence"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for edit in edits:
        writer.writerow(edit)

"""
This script will algorithmically generate genomic edits to be applied to a reference genome in order to create edited
genomes for validation testing. The edits will be "recorded" in a CSV file before they are applied. The script will
be able to randomly generate edits within strict parameters specified by the user. Fields in the CSV that need to
generated include:

1) ID - Arbitrary identifying name

2) Edit type - Substitution, Deletion, Insertion, (consider ability to create SNPs for future use)

3) Start - Position in genome where edit sequence starts

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

# Number of edits to generate (user-defined)
num_edits = 50

# Input FASTA file (provide the actual path, eventually this will be a user-defined argument)
fasta_file = "/home/casper/PycharmProjects/fantastic-lamp/validation testing/reference genome/lambda_phage.fasta"


# Function to load the reference genome from a FASTA file
def load_reference_genome(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return record.seq, len(record.seq)


# Define edit types and their weights
edit_types = ["Substitution", "Deletion", "Insertion"]
edit_weights = [1, 1, 1]  # Equal weights for testing
bps = ["A", "T", "G", "C"]  # base pairs for random edit creation


def generate_substitution(start_pos, stop_pos, length):
    """
    This function generates an edit of the substitution type. The reference sequence in the edit range will be
    replaced by a different sequence that for testing purposes will be a sequence of "N"s. Will take
    edit parameters as arguments and use these to create the edit.
    """
    es = ""  # edit sequence
    bp = random.choices(bps, k=length)
    es = es.join(bp)
    return es


def generate_deletion(start_pos, stop_pos, length):
    """
    This function generates a deletion edit. The reference sequence in the edit range will be deleted. Will take
    edit parameters as arguments and use these to create the edit.
    """
    # Determine the length of the deletion
    deletion_length = stop_pos - start_pos + 1

    # Create a string of '-' characters of the same length as the deletion
    deletion_sequence = '-' * deletion_length

    return deletion_sequence


def generate_insertion(start_pos, stop_pos, length):
    """
    This function generates an insertion edit. After the first base at the start position, a new edit sequence will
    be inserted. For testing, the sequence will consist of "N"s. Will take
    edit parameters as arguments and use these to create the edit.
    """
    es = ""  # edit sequence
    bp = random.choices(bps, k=length)
    es = es.join(bp)
    return es


# Function to generate a random edit record

def generate_edit(chromosome_length):
    edit = {
        "ID": "ID-XXXX",  # Will generate sequential IDs
        "Edit type": random.choices(edit_types, weights=edit_weights)[0],
        "Start": random.randint(100, chromosome_length - 100),
    }
    # Generate a random length within a specific range
    edit_length = random.randint(6, 60)  # Adjust the range as needed
    # Calculate the stop position based on the start and length
    edit["Stop"] = edit["Start"] + edit_length
    return edit




# Load reference genome and its length
reference_sequence, chromosome_length = load_reference_genome(fasta_file)


def generate_homology_arm_sequences(reference_sequence, edit_type, start_pos, stop_pos, edit_sequence):
    """
    This function generates the HomologyArmRefSeq and HomologyArmEditSeq fields based on the edit type and positions.
    """
    homology_arm_ref = 0
    homology_arm_edit = 0
    # if edit_type == "Deletion":
    #     field_len = 138
    # else:
    #     field_len = 138 - len(edit_sequence)
    field_len = 138
    flank_length = field_len // 2  # You can adjust this value based on your requirements

    if edit_type == "Insertion":
        # Extract the left and right flanking sequences from the reference genome
        left_flank = str(reference_sequence[start_pos - flank_length - 1: start_pos - 1])  # Adjust for 0-based indexing
        right_flank = str(reference_sequence[stop_pos: stop_pos + flank_length])
        # For insertions, the edit sequence is inserted between the flanking sequences
        homology_arm_ref = left_flank + '-' * len(edit_sequence) + right_flank
        homology_arm_edit = left_flank + edit_sequence + right_flank

    elif edit_type == "Substitution":
        # Extract the left and right flanking sequences from the reference genome
        left_flank = str(reference_sequence[start_pos - flank_length - 1: start_pos - 1])  # Adjust for 0-based indexing
        right_flank = str(reference_sequence[stop_pos: stop_pos + flank_length])
        # For substitutions and the edit sequence replaces the reference sequence
        homology_arm_ref = left_flank + reference_sequence[start_pos - 1: stop_pos] + right_flank
        homology_arm_edit = left_flank + edit_sequence + right_flank

    elif edit_type == "Deletion":
        # For deletions, extend the flanking sequences to meet the required lengths
        # Extract the left and right flanking sequences from the reference genome
        left_flank = str(reference_sequence[start_pos - flank_length - 1 - len(edit_sequence): start_pos - 1])
        # Adjusted for 0-based indexing
        if field_len % 2 == 0:  # This ensures that the length of the fields are consistent.
            right_flank = str(reference_sequence[stop_pos + len(edit_sequence): stop_pos + flank_length])
        else:
            right_flank = str(reference_sequence[stop_pos + len(edit_sequence): stop_pos + flank_length + 1])

        homology_arm_ref = left_flank + reference_sequence[start_pos:stop_pos] + right_flank
        homology_arm_edit = left_flank + right_flank

    return homology_arm_ref, homology_arm_edit


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

    homology_arm_ref, homology_arm_edit = generate_homology_arm_sequences(reference_sequence, edit["Edit type"],
                                                                          start, stop, edit_sequence)

    edit["ID"] = _ + 1
    edit["Ref sequence"] = str(reference_sequence[start - 1: stop])  # Adjust for 0-based indexing
    if edit["Edit type"] == "Deletion":
        edit["Edit sequence"] = ""
    else:
        edit["Edit sequence"] = edit_sequence
    edit["HomologyArmRefSeq"] = homology_arm_ref  # THE REF SEQ IS NOT CONTAINED WITHIN
    edit["HomologyArmEditSeq"] = homology_arm_edit
    # edit["Edit length"] = len(edit_sequence)
    edits.append(edit)

# Write edits to a CSV file
with open("genomic_edits.csv", "w", newline="") as csvfile:
    fieldnames = ["ID", "Edit type", "Start", "Stop", "Ref sequence", "Edit sequence", "HomologyArmRefSeq",
                  "HomologyArmEditSeq"]
    # fieldnames = ["ID", "Edit type", "Start", "Stop", "Ref sequence", "Edit sequence","Edit length",
    # "HomologyArmRefSeq", "HomologyArmEditSeq"] for debugging
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for edit in edits:
        writer.writerow(edit)

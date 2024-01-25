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

c) Types of edits - Substitution, Deletion, or Insertion. Frequency for each can be a weight selected. For testing, we
    can make all of them equally likely.

d) Number of edits - the user can specify how many edits should be generated.
"""

import csv
import random
from Bio import SeqIO

# Number of edits to generate (user-defined)
num_edits = 100

# Input FASTA file (provide the actual path, eventually this will be a user-defined argument)
# fasta_file = "/home/casper/PycharmProjects/fantastic-lamp/validation testing/reference genome/lambda_phage.fasta"
fasta_file = "/home/casper/PycharmProjects/fantastic-lamp/validation testing/reference genome/lambda_phage.fasta"


# Function to load the reference genome from a FASTA file
def load_reference_genome(fasta_file):
    """
    This function loads the reference genome from a FASTA file and returns the sequence as a string, and its length.
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(len(record.seq))
        return record.seq, len(record.seq)


# Define edit types and their weights
edit_types = ["Substitution", "Deletion", "Insertion"]

edit_weights = [1, 1, 1]  # These weights can be adjusted to change the frequency of each edit type

bps = ["A", "T", "G", "C"]


# base pairs for random edit creation, can be adjusted to include Ns or other characters for testing


def generate_substitution(start_pos, stop_pos, length):
    """
    This function generates an edit of the substitution type. The reference sequence in the edit range will be
    replaced by a different sequence of the same length.
    """
    es = ""  # edit sequence
    bp = random.choices(bps, k=length)  # Creates a random sequence of base pairs to use as the edit sequence
    es = es.join(bp)
    return es


def generate_deletion(start_pos, stop_pos, length):
    """
    This function generates a deletion edit. The reference sequence in the edit range will be deleted. Will take
    edit parameters as arguments and use these to create the edit. The output of this function is NOT currently used.
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
    bp = random.choices(bps, k=length)  # Creates a random sequence of base pairs to use as the edit sequence
    es = es.join(bp)
    return es


# Function to generate a random edit record
def generate_edit(chromosome_length, start_pos):
    """
    This function generates a random edit record based on the edit parameters.
    For testing, edits will be created 200bp apart.
    """

    edit_min = 20  # Minimum length of edit
    edit_max = 25  # Maximum length of edit
    margin = 200  # This is how close to the start or end of the chromosome an edit can be. For obvious reasons,
    # the edit cannot be closer to the start or end of a chromosome than the length of the edit itself.

    edit = {
        "ID": "ID-XXXX",  # Sequential IDs are added in a function lower down.
        "Edit type": random.choices(edit_types, weights=edit_weights)[0],
        "Start": start_pos
    }

    # Generate a random length within a specific range
    edit_length = random.randint(edit_min, edit_max)  # Adjust the range as needed
    # Calculate the stop position based on the start and length
    edit["Stop"] = edit["Start"] + edit_length
    start_pos += 200 + edit_length
    return edit


# Load reference genome and its length
reference_sequence, chromosome_length = load_reference_genome(fasta_file)


def generate_homology_arm_sequences(reference_sequence, edit_type, start_pos, stop_pos, edit_sequence):
    """
    This function generates the HomologyArmRefSeq and HomologyArmEditSeq fields based on the edit type and positions.
    For details on how each edit is created, see the docstrings in the if statements below.
    """
    homology_arm_ref = 0
    homology_arm_edit = 0

    field_len = 138
    # Change this to adjust the length of the homology arms. In some cases they are longer than this
    # to account for the edit, but this is the minimum. The length of the homology can vary by 1bp.

    flank_length = field_len // 2  # You can adjust this value based on your requirements

    if edit_type == "Insertion":
        """
        homology_arm_ref: flank - flank (138bp)
        homology_arm_edit: flank - inserted seq - flank (138 + inserted seq length)
        """
        edit_pos = len(edit_sequence) // 2

        homology_arm_ref = str(
            reference_sequence[start_pos - flank_length + edit_pos: start_pos + flank_length + edit_pos])

        start_pos = start_pos + edit_pos

        homology_arm_edit = (str(reference_sequence[start_pos - flank_length: start_pos]) + edit_sequence +
                             str(reference_sequence[start_pos: start_pos + flank_length]))

    elif edit_type == "Substitution":
        """
        homology_arm_ref: flank - ref seq - flank (138bp)
        homology_arm_edit: flank - edit seq (replaces ref seq) - flank (138bp)
        """
        edit_pos = len(edit_sequence) // 2
        homology_arm_ref = str(
            reference_sequence[start_pos - flank_length + edit_pos: start_pos + flank_length + edit_pos])

        homology_arm_edit = (str(reference_sequence[start_pos - flank_length + edit_pos: start_pos - 1])
                             + edit_sequence + str(reference_sequence[stop_pos: stop_pos + flank_length - edit_pos]))

    elif edit_type == "Deletion":
        """
        homology_arm_ref: flank - ref seq - flank (138 + deletion length)
        homology_arm_edit: flank - (ref seq removed) - flank (138bp)
        """
        # If you want to extend the ref arm for deletions, shift the start and stop by the edit_pos.
        # edit_pos = len(edit_sequence) // 2
        homology_arm_ref = str(reference_sequence[start_pos - flank_length:
                                                  stop_pos + flank_length])

        homology_arm_edit = (str(reference_sequence[start_pos - flank_length: start_pos - 1]) +
                             str(reference_sequence[stop_pos: stop_pos + flank_length]))

    print(
        f"Edit type: {edit_type}, homology_arm_ref: {len(homology_arm_ref)}, "
        f"homology_arm_edit: {len(homology_arm_edit)}, start: {start_pos}, stop: {stop_pos}")
    return homology_arm_ref, homology_arm_edit


# Generate and store edit records

edits = []
start_pos = 500
for x in range(num_edits+1):
    edit = generate_edit(chromosome_length, start_pos)
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

    edit["ID"] = x
    edit["Ref sequence"] = str(reference_sequence[start - 1: stop])  # Adjusted for 0-based indexing
    if edit["Edit type"] == "Deletion":
        edit["Edit sequence"] = ""
    else:
        edit["Edit sequence"] = edit_sequence
    edit["HomologyArmRefSeq"] = homology_arm_ref
    edit["HomologyArmEditSeq"] = homology_arm_edit
    # edit["Edit length"] = len(edit_sequence)
    # for testing
    edits.append(edit)
    start_pos = stop + 200

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


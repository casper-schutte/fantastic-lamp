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
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Number of edits to generate (user-defined)
num_edits = int(sys.argv[1])

# Input FASTA file (provide the actual path, eventually this will be a user-defined argument)
fasta_file = sys.argv[2]


def load_reference_genome(fasta_file):
    """
    This function loads the reference genome from a FASTA file and saves each chromosome separately.
    It returns a list of sequences and their lengths.
    """
    sequences = []
    lengths = []
    chromosomes = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(record.seq)
        lengths.append(len(record.seq))
        chromosome_id = record.id
        chromosomes[chromosome_id] = record.seq
    return sequences, lengths, chromosomes


# Load reference genome and its length
reference_sequence, chromosome_length, chrom_dict = load_reference_genome(fasta_file)
chrome_names = list(chrom_dict.keys())
# print(f"chromosome names: {chrome_names}")


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


start_positions = {name: [] for name in chrome_names}
# print(start_positions)


def find_valid_start_position(chromosome_length, existing_positions, edit_length, margin):
    """
    Find a valid start position that is at least margin away from existing positions.
    """
    valid_positions = [pos for pos in range(200, chromosome_length - 200 - edit_length) if
                       all(abs(pos - p) >= margin for p in existing_positions)]
    return valid_positions


def generate_edit(chromosome_length, chromosome_start_positions):
    """
    This function generates a random edit record based on the edit parameters.
    Edits will be created at least 200bp apart.
    """
    chrm = random.choice(chrome_names)
    edit_min = 20  # Minimum length of edit
    edit_max = 25  # Maximum length of edit
    margin = 200  # This is how close to the start or end of the chromosome an edit can be.

    if start_positions[chrm] is None:
        start_positions[chrm] = []

    while True:
        # Generate a random length within a specific range
        edit_length = random.randint(edit_min, edit_max)  # Adjust the range as needed

        # Find all valid start positions that satisfy the distance constraints
        valid_positions = find_valid_start_position(len(chrom_dict.get(chrm)), start_positions[chrm], edit_length,
                                                    margin)

        if not valid_positions:
            # If no valid positions are found, increase the margin and try again
            margin += 50
            continue

        new_start_pos = random.choice(valid_positions)
        break

    stop_pos = new_start_pos + edit_length

    edit = {
        "ID": "ID-XXXX",  # Sequential IDs are added in a function lower down.
        "Edit type": random.choices(edit_types, weights=edit_weights)[0],
        "Start": new_start_pos,
        "Chromosome": chrm  # Use the selected chromosome
    }

    # Calculate the stop position based on the start and length
    edit["Stop"] = stop_pos
    updated_start_pos = edit["Stop"]  # Update the starting position
    start_positions[edit["Chromosome"]].append(updated_start_pos)
    return edit, updated_start_pos


def generate_homology_arm_sequences(reference_sequence, edit_type, start_pos, stop_pos, edit_sequence, edit_id):
    """
    This function generates the HomologyArmRefSeq and HomologyArmEditSeq fields based on the edit type and positions.
    For details on how each edit is created, see the docstrings in the if statements below.
    """
    reference_sequence = str(reference_sequence)
    homology_arm_ref = 0
    homology_arm_edit = 0
    full_edited_genome = chrom_dict.copy()

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

        edited_genome = (str(reference_sequence[: start_pos]) + edit_sequence +
                         str(reference_sequence[start_pos:]))

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

        edited_genome = (str(reference_sequence[: start_pos - 1]) + edit_sequence +
                         str(reference_sequence[stop_pos:]))

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

        edited_genome = (str(reference_sequence[: start_pos - 1]) + str(reference_sequence[stop_pos:]))

    full_edited_genome[edit["Chromosome"]] = edited_genome

    other_chromosomes_records = [SeqRecord(Seq(chrom_sequence), id=chromosome_id, description="")
                                 for chromosome_id, chrom_sequence in full_edited_genome.items()]

    # full_edited_genome = chrom_dict
    with open(f"genome{edit_id}.fasta", "w") as fasta_file:
        if edit_id == 0:
            negative_control = SeqRecord(Seq(reference_sequence), id="negative_control", description="")
            SeqIO.write(negative_control, fasta_file, "fasta")
        else:
            SeqIO.write(other_chromosomes_records, fasta_file, "fasta")

    return homology_arm_ref, homology_arm_edit


#  Generate and store edit records
edits = []

for x in range(num_edits):
    edit, updated_start_pos = generate_edit(chromosome_length, start_positions)
    edit["ID"] = x
    start_positions[edit["Chromosome"]] = start_positions[edit["Chromosome"]].append(updated_start_pos)
    # print(start_positions)
    start = edit["Start"]
    stop = edit["Stop"]

    if edit["Edit type"] == "Substitution":
        edit_sequence = generate_substitution(start, stop, stop - start + 1)
    elif edit["Edit type"] == "Deletion":
        edit_sequence = generate_deletion(start, stop, stop - start + 1)
    elif edit["Edit type"] == "Insertion":
        edit_sequence = generate_insertion(start, stop, stop - start + 1)

    homology_arm_ref, homology_arm_edit = generate_homology_arm_sequences(chrom_dict.get(edit["Chromosome"]),
                                                                          edit["Edit type"],
                                                                          start, stop, edit_sequence, x)

    edit["Ref sequence"] = str(chrom_dict.get(edit["Chromosome"])[start - 1: stop])  # Adjusted for 0-based indexing
    if edit["Edit type"] == "Deletion":
        edit["Edit sequence"] = ""
    else:
        edit["Edit sequence"] = edit_sequence
    edit["HomologyArmRefSeq"] = homology_arm_ref
    edit["HomologyArmEditSeq"] = homology_arm_edit

    edits.append(edit)

# Write edits to a CSV file
with open("genomic_edits.csv", "w", newline="") as csvfile:
    fieldnames = ["ID", "Edit type", "Chromosome", "Start", "Stop", "Ref sequence", "Edit sequence",
                  "HomologyArmRefSeq", "HomologyArmEditSeq"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for edit in edits:
        writer.writerow(edit)

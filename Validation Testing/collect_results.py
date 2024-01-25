import csv
import sys
"""
This script will collect the relevant results from the .tsv files used for validation testing and write them to a new
file. This will allow the results to be easily compared to the expected results. The relevant results are:
- The name of the edit (from the file name)
- The edit homology arm with the highest fractional coverage
- The fractional coverage for the above edit homology arm.
"""
threshold = 0.5


def read_tsv_file(file_name):
    """
    Reads a .tsv file and returns the relevant information.
    """
    try:
        with open(file_name, newline='') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            tsv_list = list(reader)[1:6]

            # Filter non-zero rows using list comprehension
            non_zero = [row for row in tsv_list if float(row[1]) != 0]

            # Extract genome number from the file name
            genome_number = int(file_name.split("genome")[1].split(".")[0])

            # Check condition based on the extracted genome number
            correct = non_zero and genome_number == int(non_zero[0][0][13:-1]) and float(non_zero[0][1]) > threshold

    except FileNotFoundError:
        print(f"File {file_name} not found.")
        return None
    print(file_name)
    # return file_name, correct, non_zero
    if non_zero:
        return file_name, correct, non_zero[0]
    else:
        return file_name, correct


# Example usage
# for i in range(0, 101):
#     print(read_tsv_file(f"genome{i}.tsv"))

# Create a CSV file to write the results
output_file = f"{sys.argv[1]}.csv"

with open(output_file, mode='w', newline='') as csv_file:
    fieldnames = ['File Name', 'Correct', 'Edit Homology Arm', 'Fractional Coverage']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    # Write the header
    writer.writeheader()

    # Write the data for each genome
    print(f"sys-arg: {sys.argv[2]}")
    for i in range(0, int(sys.argv[2])):
        result = read_tsv_file(f"genome{i}.tsv")
        if len(result) == 2:
            writer.writerow({
                'File Name': result[0][:-4],
                'Correct': "False",
                'Edit Homology Arm': "",
                'Fractional Coverage': ""
            })
        else:
            fractional_coverage = round(float(result[2][1]), 6)
            writer.writerow({
                'File Name': result[0][:-4],
                'Correct': result[1],
                'Edit Homology Arm': result[2][0],
                'Fractional Coverage': fractional_coverage
            })

print(f"Results written to {output_file}")

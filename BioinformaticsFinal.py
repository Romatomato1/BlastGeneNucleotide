import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Blast import NCBIWWW, NCBIXML
from tabulate import tabulate


def read_file_contents(file_path):
    """
    Read the contents of a file and return as a string.
    """
    try:
        with open(file_path, 'r') as file:
            file_contents = file.read()
        return file_contents
    except IOError:
        print(f"Error: Unable to read file {file_path}")
        return ""

fasta_file1 = r"C:\Users\roman\PycharmProjects\pythonProject3\NBPF10.fasta"
fasta_file2 = r"C:\Users\roman\PycharmProjects\pythonProject3\NBPF11.fasta"

# Read the contents of the FASTA files
gene1_seq = read_file_contents(fasta_file1)
gene2_seq = read_file_contents(fasta_file2)

# Calculate the differences in nucleotides between the sequences
differences = sum(a != b for a, b in zip(gene1_seq, gene2_seq))

# Display the number of differences
print("Number of differences without BLAST:", differences)

# Create a bar plot to visualize the differences in nucleotides without BLAST
nucleotides = ['A', 'T', 'C', 'G']
gene1_counts = [gene1_seq.count(nucleotide) for nucleotide in nucleotides]
gene2_counts = [gene2_seq.count(nucleotide) for nucleotide in nucleotides]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
sns.barplot(x=nucleotides, y=gene1_counts, ax=ax1)
ax1.set_xlabel("Nucleotides")
ax1.set_ylabel("Count")
ax1.set_title("NBPF10 Nucleotide Counts (without BLAST)")

sns.barplot(x=nucleotides, y=gene2_counts, ax=ax2)
ax2.set_xlabel("Nucleotides")
ax2.set_ylabel("Count")
ax2.set_title("NBPF11 Nucleotide Counts (without BLAST)")

plt.tight_layout()
plt.show()

# Perform BLAST comparison for gene 1
result_handle_gene1 = NCBIWWW.qblast("blastn", "nt", gene1_seq)

# Parse the BLAST result for gene 1
blast_record_gene1 = NCBIXML.read(result_handle_gene1)

# Get the alignment results for gene 1
alignments_gene1 = blast_record_gene1.alignments

# Perform BLAST comparison for gene 2
result_handle_gene2 = NCBIWWW.qblast("blastn", "nt", gene2_seq)

# Parse the BLAST result for gene 2
blast_record_gene2 = NCBIXML.read(result_handle_gene2)

# Get the alignment results for gene 2
alignments_gene2 = blast_record_gene2.alignments

# Create a table to store the alignment results
table_data = []

# Display the alignment results for gene 1
for alignment_gene1 in alignments_gene1:
    alignment_info = []
    alignment_info.append(alignment_gene1.title)
    alignment_info.append(alignment_gene1.length)
    alignment_info.append(alignment_gene1.hsps[0].expect)
    table_data.append(alignment_info)

# Display the alignment results for gene 2
for alignment_gene2 in alignments_gene2:
    alignment_info = []
    alignment_info.append(alignment_gene2.title)
    alignment_info.append(alignment_gene2.length)
    alignment_info.append(alignment_gene2.hsps[0].expect)
    table_data.append(alignment_info)

# Print the table of alignment results
headers = ["Title", "Length", "Expect"]
print(tabulate(table_data, headers=headers))

# Calculate the differences in nucleotides between the sequences
differences = sum(a != b for a, b in zip(gene1_seq, gene2_seq))

# Display the number of differences
print("Number of differences with BLAST:", differences)

# Create a bar plot to visualize the differences in nucleotides with BLAST
gene1_counts = [gene1_seq.count(nucleotide) for nucleotide in nucleotides]
gene2_counts = [gene2_seq.count(nucleotide) for nucleotide in nucleotides]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
sns.barplot(x=nucleotides, y=gene1_counts, ax=ax1)
ax1.set_xlabel("Nucleotides")
ax1.set_ylabel("Count")
ax1.set_title("NBPF10 Nucleotide Counts (with BLAST)")

sns.barplot(x=nucleotides, y=gene2_counts, ax=ax2)
ax2.set_xlabel("Nucleotides")
ax2.set_ylabel("Count")
ax2.set_title("NBPF11 Nucleotide Counts (with BLAST)")

plt.tight_layout()
plt.show()

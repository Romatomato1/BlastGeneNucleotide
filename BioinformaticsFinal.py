import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Blast import NCBIWWW, NCBIXML

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

# Display the alignment results in a table
headers = ["Sequence ID", "Length", "E-value"]
print(tabulate(table_data, headers=headers))

# Plot a bar graph of the E-values
x = range(len(table_data))
evalues = [alignment[2] for alignment in table_data]

fig, ax = plt.subplots(figsize=(12, 6))
sns.barplot(x=x, y=evalues)
ax.set_xlabel("Alignment")
ax.set_ylabel("E-value")
ax.set_title("Alignment E-values")
ax.set_xticks(x)
ax.set_xticklabels([alignment[0] for alignment in table_data], rotation=45, ha='right', fontsize=8)
plt.tight_layout()
plt.show()

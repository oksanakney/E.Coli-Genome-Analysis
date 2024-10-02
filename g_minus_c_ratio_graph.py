import matplotlib.pyplot as plt
from Bio import SeqIO


def calculate_gc_content(sequence):
    """Calculate the G-C content of a given sequence."""
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_count = len(sequence)
    if total_count > 0:
        return (g_count + c_count) / total_count
    else:
        return 0.0


def gc_content_over_windows(fasta_file, window_size):
    """Calculate G-C content over sliding windows."""
    gc_ratios = []
    positions = []

    # Read the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        genome_length = len(sequence)

        # Calculate G-C content for each window
        for start in range(0, genome_length, window_size):
            end = min(start + window_size, genome_length)  # Handle end of genome
            window_sequence = sequence[start:end]
            gc_ratio = calculate_gc_content(window_sequence)
            gc_ratios.append(gc_ratio)
            positions.append(start)

    return positions, gc_ratios


# Specify your FASTA file and window size
fasta_file_path = "E-Coli-Genome.fna"  # Update this path as needed
window_size = 500000  # Size of the sliding window in nucleotides

# Calculate G-C content over the genome
positions, gc_ratios = gc_content_over_windows(fasta_file_path, window_size)

# Plotting the G-C content
plt.figure(figsize=(12, 6))
plt.plot(positions, gc_ratios, label='G-C Ratio', color='blue')
plt.title('G-C Ratio along the E. Coli Genome')
plt.xlabel('Position in Genome (Nucleotides)')
plt.ylabel('G-C Ratio')
plt.axhline(0.5, color='red', linestyle='--', label='G-C Ratio = 0.5')
plt.legend()
plt.grid()
# Save the figure
plt.savefig('gc_ratio_plot.png', format='png')  # Save as PNG
# You can also specify other formats by changing the filename and format parameter (e.g., 'gc_ratio_plot.pdf' for PDF)
plt.show()

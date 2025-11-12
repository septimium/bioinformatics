# 1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
# 2. Implement a software application that detects repetitions (between 6b and 10b) in this DNA sequence.
# 3. Plot the frequencies of the repetitions found.
import matplotlib.pyplot as plt

def find_repetitions(dna_sequence, min_length=6, max_length=10):
    all_repetitions = {}
    sequence_upper = dna_sequence.upper()
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence_upper) - length + 1):
            pattern = sequence_upper[i:i + length]
            
            if all(nucleotide in 'ATGC' for nucleotide in pattern):
                if pattern not in all_repetitions:
                    all_repetitions[pattern] = []
                all_repetitions[pattern].append(i)
    
    repetitions = {pattern: positions for pattern, positions in all_repetitions.items() 
                   if len(positions) > 1}
    
    return repetitions


def plot_repetition_frequencies(repetitions):
    sorted_patterns = sorted(repetitions.items(), key=lambda x: len(x[1]), reverse=True)
    
    top_patterns = sorted_patterns[:20]
    
    pattern_labels = [f"{pattern}" for pattern, _ in top_patterns]
    frequencies = [len(positions) for _, positions in top_patterns]
    
    plt.figure(figsize=(12, 6))
    plt.bar(range(len(frequencies)), frequencies, color='steelblue', edgecolor='black')
    plt.xlabel('Pattern', fontsize=12)
    plt.ylabel('Frequency (Number of Occurrences)', fontsize=12)
    plt.title('Repetition Frequencies', fontsize=14, fontweight='bold')
    plt.xticks(range(len(pattern_labels)), pattern_labels, rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    plt.savefig('./repetition_frequencies.png', dpi=300, bbox_inches='tight')
    print("\nFrequency plot saved as 'repetition_frequencies.png'")
    plt.show()


def main():

    dna_sequence = """
    AGTTGAAGTTGTGGAATATCTTGTTTATGCAACATGATTGTATCATGGTTTTGGCTCAAACATGCAATTA
    TTTAGTTAACAGTACATGATATTCTTTTGTTAATGTGATAACCCCTGAGGTTCCTTTTTTACAGTGCAGG
    AAGAGGCGAAGGACGGGGGAAGATATTTAAAGAGCTTCAGACACGGATTTTGCACTTCCCTCCACAATGT
    CTGCACAGACTCTGCTGCACCTCTTGGCTCTGGTGGCCTTTTTCTTTGCAGGATCCGATGCACAACTTCA
    TGAATGTGGCTTAGCTTCTCCCAACTTCAAAATAGTTGGAGGTCAGGACGGCTCACCTGGAAGTTGGCCC
    TGGCAGGTGAAGCTTTATGGCCCGTTCGGGTGCGCAGGCTCCCTGATCAACAAAAATTGGGTTCTGACTG
    CAGCTCACTGCGTCTCTAGAACAACTCCGTCCATGTGGACGGTGGTTCTGGGGCAGCAGAATCTGAGCGA
    CACACAGATGAGAACCGGAGTGAGAAGGAACGTTGGAAGAATCATCGTACACCCCAAGTTCAATCCCTCC
    CCCATCGACAACGACATTGCGTTGCTCAGAAAAAACAACAAAATAGTAATTAATATTTAAGGAAAAACAT
    TTTAAATGTTTGATTTTAATATTTACTTGTTTGCTCTTCTGTTATT
    """

    dna_sequence = "".join(char for char in dna_sequence if char.upper() in 'ATGC')
    
    print(f"\nSequence length: {len(dna_sequence)} nucleotides")
    print("\nSearching for repetitions...")
    
    repetitions = find_repetitions(dna_sequence)
    
    print(f"\nFound {len(repetitions)} repeated patterns\n")

    sorted_patterns = sorted(repetitions.items(), 
                           key=lambda x: (len(x[1]), len(x[0])), 
                           reverse=True)
    for pattern, positions in sorted_patterns:
        print(f"{pattern} ({len(pattern)}bp): {len(positions)} occurrences at positions {positions}")
    
    print("\nGenerating frequency plot...")
    plot_repetition_frequencies(repetitions)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
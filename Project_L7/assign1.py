# 1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
# 2. Implement a software application that detects repetitions (between 6b and 10b) in this DNA sequence.

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
    print("\nDone!")

if __name__ == "__main__":
    main()
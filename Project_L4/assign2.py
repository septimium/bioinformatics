# Download from NCBI the FASTA files containing the COVID19 genome and the INFLUENZA genome. 
# Use the artificial intelligence to compare the codon frequencies between the two. 
# a. Make a chart that shows the top 10 most frequent codons for COVID 19.
# b. Make a chart that shows the top 10 most frequent codons for Influenza.
# c. Compare the two results and show the most frequent codons between the two genomes.
# d. Show in the output of the console the top 3 aminoacids for each genome.

import collections
import matplotlib.pyplot as plt

CODON_TABLE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}
# -------------------------

def read_fasta(filename: str) -> str:
    sequence = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence.append(line.strip())
        
        full_sequence = "".join(sequence).upper()
        if not full_sequence:
            print(f"Warning: No sequence data found in {filename}.")
        return full_sequence
    
    except FileNotFoundError:
        print(f"Error: File not found at {filename}")
        print("Please download the required FASTA files and save them")
        print("in the same directory as this script.")
        exit(1)
    except Exception as e:
        print(f"An error occurred reading {filename}: {e}")
        exit(1)

def count_codons(genome_sequence: str) -> collections.Counter:
    rna_sequence = genome_sequence.replace('T', 'U')
    codon_counts = collections.Counter()
    
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        
        if len(codon) == 3:
            codon_counts[codon] += 1
            
    return codon_counts

def calculate_amino_acid_counts(codon_counts: collections.Counter) -> collections.Counter:
    aa_counts = collections.Counter()
    for codon, count in codon_counts.items():
        amino_acid = CODON_TABLE.get(codon)
        if amino_acid and amino_acid != 'Stop':
            aa_counts[amino_acid] += count
            
    return aa_counts

def plot_top_codons(top_10_list: list, title: str, filename: str):
    codons, counts = zip(*top_10_list)
    
    plt.figure(figsize=(12, 7))
    plt.bar(codons, counts, color='blue')
    plt.title(title, fontsize=16)
    plt.xlabel("Codon", fontsize=12)
    plt.ylabel("Frequency (Count)", fontsize=12)
    
    for i, count in enumerate(counts):
        plt.text(i, count + 5, str(count), ha='center')
        
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Chart saved as '{filename}'")

def main():
    COVID_FILE = 'covid.fasta'
    FLU_FILE = 'influenza.fasta'

    print(f"Processing {COVID_FILE}...")
    covid_seq = read_fasta(COVID_FILE)
    covid_codons = count_codons(covid_seq)
    covid_top_10 = covid_codons.most_common(10)

    print(f"Processing {FLU_FILE}...")
    flu_seq = read_fasta(FLU_FILE)
    flu_codons = count_codons(flu_seq)
    flu_top_10 = flu_codons.most_common(10)

    print("\na) Generating chart for COVID-19...")
    plot_top_codons(
        covid_top_10,
        "Top 10 Most Frequent Codons in COVID-19 (NC_045512.2)",
        "covid_top10.png"
    )
    print("\nb) Generating chart for Influenza...")
    plot_top_codons(
        flu_top_10,
        "Top 10 Most Frequent Codons in Influenza A (H1N1 2009)",
        "influenza_top10.png"
    )
    print("\n" + ("-" * 30))
    print("c) Comparison of Top 10 Codons")
    print("-" * 30)
    
    covid_top_set = {codon for codon, count in covid_top_10}
    flu_top_set = {codon for codon, count in flu_top_10}

    print(f"COVID-19 Top 10:     {sorted(list(covid_top_set))}")
    print(f"Influenza Top 10:    {sorted(list(flu_top_set))}")
    
    common_codons = covid_top_set.intersection(flu_top_set)
    
    if common_codons:
        print(f"\nCommon in both Top 10s: {sorted(list(common_codons))}")
    else:
        print("\nCommon in both Top 10s: None")

    print("\n" + ("-" * 30))
    print("d) Top 3 Most Frequent Amino Acids")
    print("-" * 30)

    covid_aa_counts = calculate_amino_acid_counts(covid_codons)
    flu_aa_counts = calculate_amino_acid_counts(flu_codons)

    print("\nCOVID-19 Top 3 Amino Acids:")
    for aa, count in covid_aa_counts.most_common(3):
        print(f"  1. {aa}: {count} occurrences")

    print("\nInfluenza Top 3 Amino Acids:")
    for aa, count in flu_aa_counts.most_common(3):
        print(f"  1. {aa}: {count} occurrences")

if __name__ == "__main__":
    main()
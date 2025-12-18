# Download 10 influenza genomes. Adapt your application from the previous assignment in order to scan each genome for possible motifs.
# For each genome make a chart that shows the signal with the most likely locations of real functional motifs.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

exon_intron_motifs = [
    "GTCATTACTA",
    "ACACAATAGA",
    "GCGAGGGGTG",
    "GGGGGGGGGG",
    "TTTTTTTTTT",
    "AATCCAAAGA",
    "AAGAACATAA",
    "AGGGTTCAGG",
    "CTATTGTCTT"
]

nucleotide_alphabet = ['A', 'C', 'G', 'T']
nucleotide_to_index = {nucleotide: idx for idx, nucleotide in enumerate(nucleotide_alphabet)}

motif_width = len(exon_intron_motifs[0])

count_matrix = np. zeros((4, motif_width), dtype=int)
for training_sequence in exon_intron_motifs:
    for position_index, base in enumerate(training_sequence):
        count_matrix[nucleotide_to_index[base], position_index] += 1

pseudocount_value = 1
weighted_count_matrix = count_matrix + pseudocount_value
relative_frequency_matrix = weighted_count_matrix / weighted_count_matrix.sum(axis=0)

background_probability = 0.25
position_weight_matrix = np.log2(relative_frequency_matrix / background_probability)


def parse_fasta_file(file_path):
  
    sequence_data = ""
    try:
        with open(file_path, 'r', encoding='utf-8') as fasta_file:
            for line in fasta_file:
                line = line.strip()
                if line.startswith(">"):
                    # Skip header lines
                    continue
                sequence_data += line.upper()
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return ""
    return sequence_data


def calculate_motif_scores(genome_sequence, pwm_matrix, window_size):
    
    score_array = []
    sequence_length = len(genome_sequence)
    
    for start_pos in range(sequence_length - window_size + 1):
        window_sequence = genome_sequence[start_pos:start_pos + window_size]
        window_score = 0.0
        
        for position, nucleotide in enumerate(window_sequence):
            if nucleotide in nucleotide_to_index: 
                window_score += pwm_matrix[nucleotide_to_index[nucleotide], position]
        
        score_array.append(window_score)
    
    return score_array


def identify_top_scoring_regions(score_list, num_peaks=5):
    
    scores_array = np.array(score_list)
    
    top_peak_indices = scores_array.argsort()[-num_peaks:][::-1]
    top_peak_scores = scores_array[top_peak_indices]
    
    return top_peak_indices, top_peak_scores


figure, subplot_axes = plt.subplots(5, 2, figsize=(16, 22))
subplot_axes = subplot_axes.flatten()

number_of_genomes = 10
successfully_processed = 0

for genome_index in range(1, number_of_genomes + 1):
    fasta_filename = f"Influenza{genome_index}.fasta"
    
    print(f"Processing {fasta_filename}...")
    
    if not os.path.exists(fasta_filename):
        print(f"  ⚠ Warning: {fasta_filename} not found, skipping...")
        print()
        
        ax = subplot_axes[genome_index - 1]
        ax.text(0.5, 0.5, f'{fasta_filename}\nNOT FOUND', 
                ha='center', va='center', fontsize=12, color='red')
        ax.set_title(f"Influenza Genome {genome_index} - Missing Data")
        ax.axis('off')
        continue
    
    genome_sequence = parse_fasta_file(fasta_filename)
    
    if not genome_sequence:
        print(f"  ⚠ Error: Could not read sequence from {fasta_filename}")
        print()
        continue
    
    motif_scores = calculate_motif_scores(genome_sequence, position_weight_matrix, motif_width)
    
    candidate_positions, candidate_scores = identify_top_scoring_regions(motif_scores, num_peaks=5)
    
    print(f"  Sequence length: {len(genome_sequence)} bp")
    print(f"  Top 5 motif candidates found at positions: {candidate_positions. tolist()}")
    print(f"  Corresponding scores: {[f'{score:.3f}' for score in candidate_scores]}")
    print()
    
    ax = subplot_axes[genome_index - 1]
    
    ax.plot(motif_scores, color='steelblue', linewidth=0.8, alpha=0.7, 
            label='Motif Score Landscape')
    
    ax.scatter(candidate_positions, candidate_scores, 
               color='crimson', s=80, zorder=5, marker='o',
               edgecolors='darkred', linewidths=1.5,
               label='Top 5 Candidate Motifs')
    
    ax.set_title(f"Influenza Genome {genome_index} - Exon-Intron Boundary Signals", 
                 fontsize=11, fontweight='bold')
    ax.set_xlabel("Genomic Position (bp)", fontsize=9)
    ax.set_ylabel("Log-Likelihood Score", fontsize=9)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='upper right', fontsize=8)
    
    mean_score = np.mean(motif_scores)
    ax.axhline(y=mean_score, color='gray', linestyle=':', linewidth=1, 
               alpha=0.5, label=f'Mean Score')
    
    successfully_processed += 1

print("=" * 80)
print(f"ANALYSIS COMPLETE: {successfully_processed}/{number_of_genomes} genomes processed successfully")
print("=" * 80)
print()


plt.tight_layout()
plt.savefig('influenza_motif_analysis.png', dpi=300, bbox_inches='tight')
print("\n Visualization saved as 'influenza_motif_analysis.png'")
plt.show()
# A very early step in splice site recognition is exon definition, a process that is as yet poorly understood. Communication between the two ends of an exon is thought to be required for this step. Computational discovery of the exon-intron border or the intron-exon border or the transcription factor binding sites (TFBS) is a challenging but important problem of bioinformatics. Implement a software application for DNA motif finding by following the steps below.
# These sequences represent the exon-intron boundary.
# 1. make the count matrix
# 2. make the weight matrix
# 3. make the relative frequencies matrix
# 4. make the Log-likelihoods Matrix
# 5. Analize sequence S by using the Log-likelihoods Matrix:
# S="CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"
# Calculate the score for each sliding window.
# Do you have signals indicating that the S sequence contains an exon-intron border?

import numpy as np
import pandas as pd

motif_sequences = [
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

motif_width = len(motif_sequences[0])
total_motifs = len(motif_sequences)

count_matrix = np.zeros((4, motif_width), dtype=int)

for sequence in motif_sequences:
    for pos, base in enumerate(sequence):
        count_matrix[nucleotide_to_index[base], pos] += 1

count_dataframe = pd.DataFrame(
    count_matrix,
    index=nucleotide_alphabet,
    columns=range(1, motif_width + 1)
)

pseudocount_value = 1
weighted_count_matrix = count_matrix + pseudocount_value

total_observations = weighted_count_matrix.sum(axis=0)
relative_frequency_matrix = weighted_count_matrix / total_observations

frequency_dataframe = pd.DataFrame(
    relative_frequency_matrix,
    index=nucleotide_alphabet,
    columns=range(1, motif_width + 1)
)

background_probability = 0.25  
log_likelihood_matrix = np.log2(relative_frequency_matrix / background_probability)

log_likelihood_dataframe = pd. DataFrame(
    log_likelihood_matrix,
    index=nucleotide_alphabet,
    columns=range(1, motif_width + 1)
)

target_sequence = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

sliding_windows = []
window_scores = []

for start_position in range(len(target_sequence) - motif_width + 1):
    current_window = target_sequence[start_position: start_position + motif_width]
    
    window_score = sum(
        log_likelihood_matrix[nucleotide_to_index[base], col_idx]
        for col_idx, base in enumerate(current_window)
    )
    
    sliding_windows.append(current_window)
    window_scores.append(window_score)

scores_dataframe = pd.DataFrame({
    "Position": range(len(sliding_windows)),
    "Window_Sequence": sliding_windows,
    "Log_Likelihood_Score": window_scores
})

print("=" * 70)
print("STEP 1: COUNT MATRIX")
print("=" * 70)
print(count_dataframe.to_string())
print()

print("=" * 70)
print("STEP 3: RELATIVE FREQUENCIES MATRIX")
print("=" * 70)
print(frequency_dataframe.to_string(float_format="%.4f"))
print()

print("=" * 70)
print("STEP 4: LOG-LIKELIHOODS MATRIX")
print("=" * 70)
print(log_likelihood_dataframe. to_string(float_format="%.4f"))
print()

print("=" * 70)
print("STEP 5: SLIDING WINDOW ANALYSIS OF SEQUENCE S")
print("=" * 70)
print(scores_dataframe.to_string(index=False, float_format="%.4f"))
print()

max_score = max(window_scores)
max_score_position = window_scores. index(max_score)
max_score_window = sliding_windows[max_score_position]

print("=" * 70)
print("ANALYSIS SUMMARY")
print("=" * 70)
print(f"Highest scoring window: '{max_score_window}' at position {max_score_position}")
print(f"Maximum log-likelihood score: {max_score:.4f}")
print()

print("INTERPRETATION:")
print("-" * 70)
if max_score > 0:
    print("YES - Signal detected!")
    print(f"The sequence S exhibits significant matches to the exon-intron")
    print(f"boundary motif, with a positive log-likelihood score of {max_score:.4f}")
    print(f"at position {max_score_position}.  This positive score indicates that")
    print(f"the subsequence '{max_score_window}' is more likely to belong to")
    print(f"the exon-intron border pattern than to random background sequence,")
    print(f"suggesting the presence of a splice site recognition signal.")
else:
    print("NO significant signal detected.")
    print(f"The maximum score ({max_score:.4f}) is not sufficiently high")
    print(f"to confidently identify an exon-intron boundary in sequence S.")
print("=" * 70)
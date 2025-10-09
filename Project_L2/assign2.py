# Find in sequence S only the dinucleotides and trinucleotides that exist, without the use of the bruteforce engine.
# In order to achieve the results, one must verify these combinations starting from the beginning of the sequence.

# ex: S = "abaa"

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

def find_existing_combinations(sequence, length):
    total = len(sequence) - length + 1
    counts = {}

    for i in range(total):
        sub = sequence[i:i + length]
        if sub in counts:
            counts[sub] += 1
        else:
            counts[sub] = 1

    percentages = {combo: (count / total) * 100 for combo, count in counts.items()}
    return percentages

di_percentages = find_existing_combinations(S, 2)

tri_percentages = find_existing_combinations(S, 3)

print("Dinucleotide Percentages from Sequence: ")
for combo, perc in sorted(di_percentages.items()):
    print(f"{combo}: {perc:.2f}%")

print("Trinucleotide Percentages from Sequence: ")
for combo, perc in sorted(tri_percentages.items()):
    print(f"{combo}: {perc:.2f}%")
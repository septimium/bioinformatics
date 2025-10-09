# Find the percentage for all the dinucleotide and trinucleotide combinations for the sequence: S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA".

# 1. Build a brute force engine to generate all dinucleotide and trinucleotide combinations.
# 2. For each combination, find out the percentage inside the S sequence.
# 3. Show the percentage for each combination in the output of your implementation.
    
S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

def generate_2combinations(letters):
    combos = []
    for l1 in letters:
        for l2 in letters:
            combos.append(l1 + l2)
    return combos

def generate_3combinations(letters):
    combos = []
    for l1 in letters:
        for l2 in letters:
            for l3 in letters: 
                combos.append(l1 + l2 + l3)
    return combos

def count_combinations_and_percentages(sequence, length):
    letters = ['A', 'C', 'G', 'T']
    if(length == 2):
        combinations = generate_2combinations(letters)
    if(length == 3):
        combinations = generate_3combinations(letters)
    total = len(sequence) - length + 1
    counts = {combo: 0 for combo in combinations}

    for combo in combinations:
        for i in range(total):
            match = True
            for j in range(length):
                if sequence[i + j] != combo[j]:
                    match = False
                    break
            if match:
                counts[combo] += 1

    percentages = {}
    for combo in combinations:
        percentages[combo] = (counts[combo] / total) * 100

    return percentages


two_letter = count_combinations_and_percentages(S, 2)
three_letter = count_combinations_and_percentages(S, 3)

print("Dinucleotide Percentages from Sequence: ")
for combo in sorted(two_letter.keys()):
    if(two_letter[combo] != 0):
        print(f"{combo}: {two_letter[combo]:.2f}%")

print("Trinucleotide Percentages from Sequence: ")
for combo in sorted(three_letter.keys()):
    if(three_letter[combo] != 0):
        print(f"{combo}: {three_letter[combo]:.2f}%")

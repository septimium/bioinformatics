# Use a DNA sequence of 50 letters in order to calculate the transition matrix that represents this sequence and store it in a JSON file.
import json
import numpy as np
from collections import defaultdict

def generate_dna_sequence(length=50):
    import random
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def calculate_transition_matrix(sequence):
    bases = ['A', 'C', 'G', 'T']
    
    transition_counts = defaultdict(lambda: defaultdict(int))
    
    for i in range(len(sequence) - 1):
        current_base = sequence[i]
        next_base = sequence[i + 1]
        transition_counts[current_base][next_base] += 1
    
    transition_matrix = {}
    for base in bases:
        total = sum(transition_counts[base].values())
        if total > 0:
            transition_matrix[base] = {
                next_base: transition_counts[base][next_base] / total
                for next_base in bases
            }
        else:
            transition_matrix[base] = {next_base: 0.25 for next_base in bases}
    
    matrix_array = np.array([
        [transition_matrix[from_base][to_base] for to_base in bases]
        for from_base in bases
    ])
    
    return transition_matrix, matrix_array

def save_to_json(data, filename="assign2.json"):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    print(f"\nTransition matrix saved to '{filename}'")

def print_transition_matrix(matrix_dict, matrix_array):
    bases = ['A', 'C', 'G', 'T']
    
    print("\nTransition Matrix (Dictionary Format):")
    print("="*60)
    for from_base in bases:
        print(f"\nFrom {from_base}:")
        for to_base in bases:
            prob = matrix_dict[from_base][to_base]
            print(f"  → {to_base}: {prob:.4f} ({prob*100:.2f}%)")
    
    print("\n" + "="*60)
    print("\nTransition Matrix (Array Format):")
    print("     ", "  ".join(bases))
    print("    " + "-"*25)
    for i, from_base in enumerate(bases):
        row_str = f"{from_base} | "
        row_str += "  ".join(f"{matrix_array[i][j]:.3f}" for j in range(4))
        print(row_str)

def main():
    
    print("DNA Sequence Transition Matrix Generator")
    print("="*60)
    
    dna_sequence = generate_dna_sequence(50)
    
    print(f"\nDNA Sequence (length: {len(dna_sequence)}):")
    print(dna_sequence)
    
    transition_dict, transition_array = calculate_transition_matrix(dna_sequence)
    
    print_transition_matrix(transition_dict, transition_array)
    
    json_data = {
        "dna_sequence": dna_sequence,
        "sequence_length": len(dna_sequence),
        "transition_matrix": transition_dict,
        "matrix_array": transition_array.tolist(),
        "bases": ['A', 'C', 'G', 'T'],
        "description": "Transition probabilities between nucleotides in DNA sequence"
    }
    
    save_to_json(json_data)
    
    print("\n" + "="*60)
    print("\nSequence Statistics:")
    bases = ['A', 'C', 'G', 'T']
    for base in bases:
        count = dna_sequence.count(base)
        percentage = (count / len(dna_sequence)) * 100
        print(f"  {base}: {count} ({percentage:.2f}%)")
    
    print("\n" + "="*60)
    print("\nMatrix Properties:")
    print(f"  Sum of each row (should be ~1.0):")
    for i, base in enumerate(bases):
        row_sum = np.sum(transition_array[i])
        print(f"    Row {base}: {row_sum:.6f}")

if __name__ == "__main__":
    main()
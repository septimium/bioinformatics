import numpy as np
from collections import defaultdict

# Known sequences
S1 = "ATCGATTCGATATCATACACGTAT"  # CpG island (+)
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"  # Non-CpG island (-)
S = "CAGGTTGGAAACGTAA"  # Test sequence

nucleotides = ['A', 'C', 'G', 'T']

def count_transitions(sequence):
    """Count transition frequencies from one nucleotide to another"""
    transitions = defaultdict(lambda: defaultdict(int))
    
    for i in range(len(sequence) - 1):
        current = sequence[i]
        next_nt = sequence[i + 1]
        transitions[current][next_nt] += 1
    
    return transitions

def calculate_probabilities(transitions):
    """Convert transition counts to probabilities"""
    probabilities = {}
    
    for current in nucleotides:
        total = sum(transitions[current][next_nt] for next_nt in nucleotides)
        probabilities[current] = {}
        
        for next_nt in nucleotides:
            if total > 0:
                probabilities[current][next_nt] = transitions[current][next_nt] / total
            else:
                probabilities[current][next_nt] = 0.25  # uniform distribution
    
    return probabilities

def create_log_likelihood_matrix(prob_plus, prob_minus):
    """Create log-likelihood matrix: log(P+/P-)"""
    log_likelihood = {}
    
    for current in nucleotides:
        log_likelihood[current] = {}
        for next_nt in nucleotides:
            p_plus = prob_plus[current][next_nt]
            p_minus = prob_minus[current][next_nt]
            
            # Add pseudocount to avoid log(0)
            p_plus = max(p_plus, 0.01)
            p_minus = max(p_minus, 0.01)
            
            log_likelihood[current][next_nt] = np.log(p_plus / p_minus)
    
    return log_likelihood

def calculate_sequence_score(sequence, log_likelihood):
    """Calculate total log-likelihood score for a sequence"""
    score = 0
    
    for i in range(len(sequence) - 1):
        current = sequence[i]
        next_nt = sequence[i + 1]
        if current in log_likelihood and next_nt in log_likelihood[current]:
            score += log_likelihood[current][next_nt]
    
    return score

# Step 1: Count transitions for CpG+ model (S1)
print("Step 1: CpG+ Model (S1)")
print(f"Sequence S1: {S1}")
transitions_plus = count_transitions(S1)
print("\nTransition counts (S1):")
for current in nucleotides:
    for next_nt in nucleotides:
        count = transitions_plus[current][next_nt]
        if count > 0:
            print(f"  {current} -> {next_nt}: {count}")

# Step 2: Count transitions for CpG- model (S2)
print("\n" + "="*50)
print("Step 2: CpG- Model (S2)")
print(f"Sequence S2: {S2}")
transitions_minus = count_transitions(S2)
print("\nTransition counts (S2):")
for current in nucleotides:
    for next_nt in nucleotides:
        count = transitions_minus[current][next_nt]
        if count > 0:
            print(f"  {current} -> {next_nt}: {count}")

# Calculate probabilities
prob_plus = calculate_probabilities(transitions_plus)
prob_minus = calculate_probabilities(transitions_minus)

# Step 3: Create log-likelihood matrix
print("\n" + "="*50)
print("Step 3: Log-Likelihood Matrix log(P+/P-)")
log_likelihood = create_log_likelihood_matrix(prob_plus, prob_minus)

print("\n    ", end="")
for next_nt in nucleotides:
    print(f"{next_nt:>8}", end="")
print()

for current in nucleotides:
    print(f"{current}  ", end="")
    for next_nt in nucleotides:
        print(f"{log_likelihood[current][next_nt]:>8.3f}", end="")
    print()

# Step 4: Classify test sequence S
print("\n" + "="*50)
print(f"Test Sequence: {S}")
score = calculate_sequence_score(S, log_likelihood)

print(f"\nLog-likelihood score: {score:.3f}")
print("\nClassification:")
if score > 0:
    print(f"✓ Sequence belongs to CpG ISLAND (score > 0)")
else:
    print(f"✗ Sequence belongs to NON-ISLAND region (score < 0)")

print(f"\nInterpretation: A positive score indicates the sequence is more likely")
print(f"to belong to a CpG island, while a negative score suggests it does not.")
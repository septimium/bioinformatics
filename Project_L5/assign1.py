import random

def get_sequence_data():
    original_sequence = (
        "GATCAATTGTCGCTTAGTTCATTACTGTTATTTTCTTTTTGTGAATATTCAATTGTTTCGA"
        "GCGGCAAAAACATATGTTTTGTTTTATCTGCTGGCTTTGTCATCTCTTAGAGATTTAGGT"
        "GTTTATATTCAAGCGTTCAGATAATAAAAAATATCTTTCAATTTTTGGCATAGCGTGTC"
        "AAATCATTTGCTTTCTCGGCAAGTGGACAGTTATCCACTAAAAATATTCAATTGGCTCA"
        "GCTTTTCATTTTCTTTTATATTTTCTTCAGCAGGGTTCAATCGTTCTTATTTTCTTTTC"
        "TGTTGATTTAGAATCTTATTTGTTTTTCTTGTTTCTGATTTAGGTAATATTTTCTCTTT"
        "ATTGATAGTAATTCTGCTTTATTGAGCTTTTCTATTTATGTTTTTTTTTCTTCTTGAGC"
        "AGTAAACAAAAAAAAATATTTTTATTTCTTTGGTTTTATTTCCTTTTTGTTTTCATTT"
        "TTTTATTTTGAATAAAACATTTTCTTATGATAAATAATTCATTTTATCAGATATTTCTT"
        "CTTTTGTTTTTTGGTATCAGTGAATTATTTTTTATTTAATAATTGCTTTTTCTTTTTAT"
        "ATTTTCTTTTTATATTATTTATATATTCAATTGGTAATAATTGTTTTCTTCTGTTTTAT"
        "TTATTGCTTGTTTTGATTATTTTTAATAATTGTTGTTTTATTTCTTTCTTATTTCTTTA"
        "TCTTATAAATCAATTTTGTTCTTTTATTCAATTTTATTTGTTATATTGTTATTTTCATT"
        "TTGCTGTTTATCAATTGCAATTGTTATCTTTTTTATGAGCAAATGATTCAATCTTTTTG"
        "TCAAAAATTTGTTTCAGTTTTAATTTTTCGATTGATGCTGTA"
    )
    return original_sequence

def sample_sequence(sequence, num_samples=2000, min_len=100, max_len=150):
    samples = []
    seq_len = len(sequence)
    for _ in range(num_samples):
        sample_length = random.randint(min_len, max_len)
        max_start_index = seq_len - sample_length
        if max_start_index < 0:
            start_index = 0
            sample_length = seq_len
        else:
            start_index = random.randint(0, max_start_index)
        sample = sequence[start_index : start_index + sample_length]
        samples.append(sample)
    return samples

def find_best_suffix_overlap(assembly, fragment, min_overlap=30):
    best_k = 0
    for k in range(min(len(assembly), len(fragment)), min_overlap - 1, -1):
        if assembly.endswith(fragment[:k]):
            best_k = k
            break  
    return best_k

def find_best_prefix_overlap(assembly, fragment, min_overlap=30):
    best_k = 0
    for k in range(min(len(assembly), len(fragment)), min_overlap - 1, -1):
        if fragment.endswith(assembly[:k]):
            best_k = k
            break 
    return best_k

def assemble_sequence(samples, min_overlap=30):
    if not samples:
        return ""
        
    remaining_samples = list(samples)
    assembly = remaining_samples.pop(0)
    print(f"Starting assembly with a {len(assembly)}bp fragment.")
    while True:
        best_overlap_len = 0
        best_match_index = -1
        
        for i, sample in enumerate(remaining_samples):
            overlap_len = find_best_suffix_overlap(assembly, sample, min_overlap)
            if overlap_len > best_overlap_len:
                best_overlap_len = overlap_len
                best_match_index = i
        
        if best_match_index != -1:
            best_sample = remaining_samples.pop(best_match_index)
            assembly += best_sample[best_overlap_len:]
        else:
            break
            
    print(f"Suffix extension complete. Assembly is now {len(assembly)}bp.")

    while True:
        best_overlap_len = 0
        best_match_index = -1

        for i, sample in enumerate(remaining_samples):
            overlap_len = find_best_prefix_overlap(assembly, sample, min_overlap)
            if overlap_len > best_overlap_len:
                best_overlap_len = overlap_len
                best_match_index = i
        
        if best_match_index != -1:
            best_sample = remaining_samples.pop(best_match_index)
            assembly = best_sample[:-best_overlap_len] + assembly
        else:
            break

    print(f"Prefix extension complete. Final assembly length: {len(assembly)}bp.")
    print(f"{len(remaining_samples)} samples were unused.")
    return assembly

def main():
    original_sequence = get_sequence_data()
    print(f"Original sequence length: {len(original_sequence)}bp")
    samples = sample_sequence(original_sequence, num_samples=2000, min_len=100, max_len=150)
    print(f"Generated {len(samples)} samples.")
    
    random.shuffle(samples) 
    
    final_assembly = assemble_sequence(samples, min_overlap=30)
    
    print("\n--- Assembly Results ---")
    print(f"Original Length:  {len(original_sequence)}")
    print(f"Assembled Length: {len(final_assembly)}")
    
    if final_assembly == original_sequence:
        print("Verification: SUCCESS! Assembly perfectly matches original.")
    elif final_assembly in original_sequence:
        print("Verification: PARTIAL. Assembly is a perfect substring of the original (likely missed ends).")
    elif original_sequence in final_assembly:
        print("Verification: PARTIAL. Original is a perfect substring of the assembly (over-assembled).")
    else:
        print("Verification: FAILED. Assembly does not match original.")


if __name__ == "__main__":
    main()


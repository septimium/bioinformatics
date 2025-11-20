# 3.⁠ ⁠Download from NCBI 3 bacterial genomes of your choosing. Try to find in these genomes possible tranposons.
# For these, one must detect possible inverted repeats without prior knowledge about their existence in the sequence.
# The inverted repeat should have 4 minimum length of 4 letters and a maximum of 6 letters.
# 4.⁠ ⁠Make a raport about the results from exercise 3 write the raport in a txt or docx file. Upload this raport on moodle.

def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):  
                sequence += line.upper()
    return sequence

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def find_inverted_repeats(sequence, min_len=4, max_len=6, max_spacer=100):
    repeats = []
    
    for repeat_len in range(min_len, max_len + 1):
        for i in range(len(sequence) - repeat_len):
            left = sequence[i:i + repeat_len]
            
            if set(left) - set('ACGT'):
                continue
            
            rc = reverse_complement(left)
            
            for j in range(i + repeat_len, min(i + repeat_len + max_spacer, len(sequence) - repeat_len + 1)):
                right = sequence[j:j + repeat_len]
                
                if right == rc:
                    spacer = j - (i + repeat_len)
                    repeats.append({
                        'left_seq': left,
                        'right_seq': right,
                        'left_pos': i,
                        'right_pos': j,
                        'length': repeat_len,
                        'spacer': spacer
                    })
    
    return repeats

def filter_repeats(repeats, max_results=50):
    """Keep only non-overlapping repeats"""
    sorted_repeats = sorted(repeats, key=lambda x: (-x['length'], -x['spacer']))
    
    filtered = []
    used_positions = set()
    
    for repeat in sorted_repeats:
        if len(filtered) >= max_results:
            break
        
        left_range = range(repeat['left_pos'], repeat['left_pos'] + repeat['length'])
        right_range = range(repeat['right_pos'], repeat['right_pos'] + repeat['length'])
        all_pos = set(left_range) | set(right_range)
        
        if not all_pos & used_positions:
            filtered.append(repeat)
            used_positions.update(all_pos)
    
    return sorted(filtered, key=lambda x: x['left_pos'])

def analyze_genome(filename):
    print(f"\nAnalyzing: {filename}")
    
    sequence = read_fasta(filename)
    print(f"  Genome length: {len(sequence):,} bp")
    
    print(f"  Searching for inverted repeats...")
    repeats = find_inverted_repeats(sequence)
    print(f"  Found {len(repeats)} total inverted repeats")
    
    filtered = filter_repeats(repeats)
    print(f"  Selected {len(filtered)} top candidates")
    
    counts = {4: 0, 5: 0, 6: 0}
    for r in filtered:
        counts[r['length']] += 1
    
    return {
        'filename': filename,
        'length': len(sequence),
        'total_repeats': len(repeats),
        'filtered_repeats': filtered,
        'counts': counts
    }

def generate_report(results):
    print("\nGenerating text report...")
    with open('transposon_report.txt', 'w') as f:
        f.write("="*80 + "\n")
        f.write("TRANSPOSON DETECTION REPORT\n")
        f.write("Analysis of Inverted Repeats\n")
        f.write("="*80 + "\n\n")
        
        f.write("METHODOLOGY:\n")
        f.write("-" * 80 + "\n")
        f.write("Transposons are detected by identifying inverted repeats (ITRs).\n")
        f.write("Inverted repeats consist of a sequence and its reverse complement.\n\n")
        f.write("Parameters:\n")
        f.write("  - Inverted repeat length: 4-6 bp\n")
        f.write("  - Maximum spacer distance: 100 bp\n\n")
        
        for i, result in enumerate(results, 1):
            f.write(f"\n{'='*80}\n")
            f.write(f"GENOME {i}: {result['filename']}\n")
            f.write(f"{'='*80}\n")
            f.write(f"Genome length: {result['length']:,} bp\n")
            f.write(f"Total inverted repeats: {result['total_repeats']}\n")
            f.write(f"Top candidates: {len(result['filtered_repeats'])}\n\n")
            
            f.write("Distribution:\n")
            for length in [4, 5, 6]:
                f.write(f"  {length} bp repeats: {result['counts'][length]}\n")
            
            f.write(f"\n--- TOP 30 INVERTED REPEATS ---\n\n")
            
            for j, repeat in enumerate(result['filtered_repeats'][:30], 1):
                f.write(f"Repeat #{j}:\n")
                f.write(f"  Left:  {repeat['left_seq']} at position {repeat['left_pos']}\n")
                f.write(f"  Right: {repeat['right_seq']} at position {repeat['right_pos']}\n")
                f.write(f"  Length: {repeat['length']} bp, Spacer: {repeat['spacer']} bp\n")
                f.write(f"  Structure: {repeat['left_seq']}--[{repeat['spacer']}bp]--{repeat['right_seq']}\n\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("CONCLUSIONS\n")
        f.write("="*80 + "\n")
        f.write("Inverted repeats were successfully detected in all genomes.\n")
        f.write("These structures are characteristic signatures of transposable elements.\n")
    
    print("  Text report saved: transposon_report.txt")

def main():
    print("="*80)
    print("TRANSPOSON DETECTION - INVERTED REPEAT ANALYSIS")
    print("="*80)
    
    print("\nEnter FASTA file paths (one per line, empty line when done):")
    fasta_files = []
    
    while True:
        filepath = input(f"Genome {len(fasta_files) + 1}: ").strip()
        if not filepath:
            break
        fasta_files.append(filepath)
    
    if not fasta_files:
        print("\nNo files provided. Exiting.")
        return
    
    print(f"\n{len(fasta_files)} genome(s) will be analyzed")
    
    results = []
    for fasta_file in fasta_files:
        try:
            result = analyze_genome(fasta_file)
            results.append(result)
        except FileNotFoundError:
            print(f"  ERROR: File not found: {fasta_file}")
        except Exception as e:
            print(f"  ERROR: {e}")
    
    if not results:
        print("\nNo genomes were successfully analyzed.")
        return
    
    generate_report(results)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print("\nReports generated:")
    print("  • transposon_report.txt")

if __name__ == "__main__":
    main()
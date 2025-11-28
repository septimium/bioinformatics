import matplotlib.pyplot as plt
import numpy as np

def find_restriction_sites(sequence, enzyme_name, recognition_seq, cut_position):
    cut_positions = []
    seq_upper = sequence.upper()
    
    start = 0
    while True:
        pos = seq_upper.find(recognition_seq, start)
        if pos == -1:
            break
        cut_pos = pos + cut_position
        cut_positions.append(cut_pos)
        start = pos + 1
    
    return sorted(cut_positions)

def digest_sequence(sequence, all_cut_positions):
    if not all_cut_positions:
        return [len(sequence)]
    
    fragments = []
    cut_positions = sorted(set(all_cut_positions))
    
    positions = [0] + cut_positions + [len(sequence)]
    
    for i in range(len(positions) - 1):
        fragment_length = positions[i + 1] - positions[i]
        if fragment_length > 0:
            fragments.append(fragment_length)
    
    return sorted(fragments, reverse=True)

def plot_gel_electrophoresis(enzyme_results, sequence_length):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8))
    
    enzyme_names = list(enzyme_results.keys())
    num_enzymes = len(enzyme_names)
    
    max_fragment = max([max(data['fragments']) if data['fragments'] else 0 
                       for data in enzyme_results.values()])
    
    ax1.set_xlim(-0.5, num_enzymes - 0.5)
    ax1.set_ylim(0, max_fragment * 1.1)
    ax1.invert_yaxis()
    ax1.set_xlabel('Enzyme', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Fragment Size (bp)', fontsize=12, fontweight='bold')
    ax1.set_title('Simulated Gel Electrophoresis', fontsize=14, fontweight='bold')
    ax1.set_xticks(range(num_enzymes))
    ax1.set_xticklabels(enzyme_names, rotation=45, ha='right')
    ax1.grid(axis='y', alpha=0.3)
    
    for i, enzyme_name in enumerate(enzyme_names):
        fragments = enzyme_results[enzyme_name]['fragments']
        for frag_size in set(fragments):
            count = fragments.count(frag_size)
            intensity = min(count / 3.0, 1.0)
            ax1.add_patch(plt.Rectangle((i - 0.3, frag_size - 20), 0.6, 40, 
                                       facecolor='blue', alpha=intensity, edgecolor='darkblue'))
    
    cut_data = [enzyme_results[name]['num_cuts'] for name in enzyme_names]
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, num_enzymes))
    bars = ax2.bar(enzyme_names, cut_data, color=colors, edgecolor='black', linewidth=1.5)
    
    ax2.set_xlabel('Enzyme', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Cuts', fontsize=12, fontweight='bold')
    ax2.set_title('Restriction Sites per Enzyme', fontsize=14, fontweight='bold')
    ax2.set_xticklabels(enzyme_names, rotation=45, ha='right')
    ax2.grid(axis='y', alpha=0.3)
    
    for bar, count in zip(bars, cut_data):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(count)}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.show()

def plot_fragment_distribution(enzyme_results):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    enzyme_names = list(enzyme_results.keys())
    
    for idx, enzyme_name in enumerate(enzyme_names):
        ax = axes[idx]
        fragments = enzyme_results[enzyme_name]['fragments']
        
        if fragments:
            ax.hist(fragments, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
            ax.set_xlabel('Fragment Size (bp)', fontsize=10)
            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_title(f'{enzyme_name}\n{len(fragments)} fragments', fontsize=11, fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No cuts', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=12)
            ax.set_title(enzyme_name, fontsize=11, fontweight='bold')
    
    axes[-1].axis('off')
    
    plt.tight_layout()
    plt.show()

def main():
    enzymes = {
        'EcoRI': {'recognition': 'GAATTC', 'cut_position': 1},
        'BamHI': {'recognition': 'GGATCC', 'cut_position': 1},
        'HindIII': {'recognition': 'AAGCTT', 'cut_position': 1},
        'TaqI': {'recognition': 'TCGA', 'cut_position': 1},
        'HaeIII': {'recognition': 'GGCC', 'cut_position': 2}
    }
    
    print("RESTRICTION ENZYME DIGESTION SIMULATOR")
    print("="*60)
    print("\nEnter the DNA sequence:")
    sequence = input().strip().upper()
    
    if not (1000 <= len(sequence) <= 3000):
        print(f"Warning: Sequence length is {len(sequence)} bp (recommended: 1000-3000 bp)")
    
    if not all(base in 'ATCG' for base in sequence):
        print("Error: Sequence contains invalid characters. Only A, T, C, G allowed.")
        return
    
    print(f"\nSequence length: {len(sequence)} bp")
    print("\nAnalyzing restriction sites...\n")
    
    enzyme_results = {}
    all_cuts = []
    
    for enzyme_name, enzyme_data in enzymes.items():
        recognition = enzyme_data['recognition']
        cut_pos = enzyme_data['cut_position']
        
        cut_positions = find_restriction_sites(sequence, enzyme_name, recognition, cut_pos)
        fragments = digest_sequence(sequence, cut_positions)
        
        enzyme_results[enzyme_name] = {
            'cut_positions': cut_positions,
            'fragments': fragments,
            'num_cuts': len(cut_positions)
        }
        
        all_cuts.extend(cut_positions)
        
        print(f"{enzyme_name}:")
        print(f"  Recognition: {recognition}")
        print(f"  Cuts: {len(cut_positions)}")
        print(f"  Cut positions: {cut_positions if len(cut_positions) <= 10 else cut_positions[:10] + ['...']}")
        print(f"  Fragments: {len(fragments)}")
        print(f"  Sizes (bp): {fragments}")
        print()
    
    combined_fragments = digest_sequence(sequence, all_cuts)
    print(f"COMBINED DIGEST:")
    print(f"  Total cuts: {len(set(all_cuts))}")
    print(f"  Fragments: {len(combined_fragments)}")
    print(f"  Sizes (bp): {combined_fragments}")
    print("\nGenerating plots...")
    
    plot_gel_electrophoresis(enzyme_results, len(sequence))
    plot_fragment_distribution(enzyme_results)

if __name__ == "__main__":
    main()
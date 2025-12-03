import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List

def calculate_cg_content(sequence: str) -> float:
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    return ((c_count + g_count) / len(sequence)) * 100

def calculate_index_of_coincidence(sequence: str) -> float:
    n = len(sequence)
    if n == 0:
        return 0.0
    
    counts = {'A': sequence.count('A'), 'C': sequence.count('C'), 
              'G': sequence.count('G'), 'T': sequence.count('T')}
    
    expected = n / 4
    chi_sq = sum((count - expected) ** 2 / expected for count in counts.values())
    
    return chi_sq

def sliding_window_analysis(sequence: str, window_size: int) -> Tuple[List[float], List[float]]:
    cg_values = []
    ic_values = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        cg_values.append(calculate_cg_content(window))
        ic_values.append(calculate_index_of_coincidence(window))
    
    return cg_values, ic_values

def process_cg_content(sequence: str, window_size: int) -> float:
    cg_values, _ = sliding_window_analysis(sequence, window_size)
    weights = np.arange(1, len(cg_values) + 1)
    return np. average(cg_values, weights=weights)

def process_index_of_coincidence(sequence: str, window_size: int) -> float:
    _, ic_values = sliding_window_analysis(sequence, window_size)
    return np.mean(ic_values)

def calculate_center_of_weight(cg_values: List[float], ic_values: List[float]) -> Tuple[float, float]:
    weights = np.arange(1, len(cg_values) + 1)
    cg_center = np.average(cg_values, weights=weights)
    ic_center = np.mean(ic_values)
    return cg_center, ic_center

def plot_pattern(cg_values: List[float], ic_values: List[float], 
                 title: str = "DNA Pattern", show_center: bool = True):
    plt.figure(figsize=(10, 6))
    
    colors = np.arange(len(cg_values))
    scatter = plt.scatter(cg_values, ic_values, c=colors, cmap='viridis', 
                         alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    plt.colorbar(scatter, label='Window Position')
    plt.xlabel('C+G Content (%)', fontsize=12)
    plt.ylabel('Index of Coincidence (χ²)', fontsize=12)
    plt. title(title, fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    if show_center:
        center_cg, center_ic = calculate_center_of_weight(cg_values, ic_values)
        plt.scatter([center_cg], [center_ic], color='red', s=300, 
                   marker='X', linewidths=3, label=f'Center ({center_cg:.2f}, {center_ic:.2f})', 
                   edgecolors='darkred', zorder=5)
        plt.legend(fontsize=10)
    
    plt. tight_layout()
    return plt.gcf()

def analyze_promoter(sequence: str, window_size: int = 30) -> dict:
    cg_values, ic_values = sliding_window_analysis(sequence, window_size)
    center = calculate_center_of_weight(cg_values, ic_values)
    
    return {
        'cg_values': cg_values,
        'ic_values': ic_values,
        'center': center,
        'avg_cg': center[0],
        'avg_ic': center[1],
        'sequence_length': len(sequence),
        'window_count': len(cg_values)
    }

if __name__ == "__main__":
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    
    window_size = 30
    
    print(f"Sequence: {S}")
    print(f"Sequence length: {len(S)}")
    print(f"Window size: {window_size}b")
    print(f"Number of windows: {len(S) - window_size + 1}\n")
    
    cg_result = process_cg_content(S, window_size)
    print(f"CG Content (Weighted Average):")
    print(f"  Expected: 29.27")
    print(f"  Computed: {cg_result:.2f}")
    
    ic_result = process_index_of_coincidence(S, window_size)
    print(f"Index of Coincidence (Average, χ²):")
    print(f"  Expected: 27.53")
    print(f"  Computed: {ic_result:.2f}")
    
    results = analyze_promoter(S, window_size)
    
    center_cg, center_ic = results['center']
    print(f"Center of Weight:")
    print(f"  CG = {center_cg:.2f}%")
    print(f"  IC = {center_ic:.2f}\n")
    
    fig1 = plot_pattern(results['cg_values'], results['ic_values'], 
                        title=f"DNA Pattern Analysis - Test Sequence ({results['window_count']} windows)")
    plt.savefig('dna_pattern.png', dpi=300, bbox_inches='tight')
    
    sequences = [
        ("Full sequence", S),
        ("5' region (1-60)", S[0:60]),
        ("Middle region (20-70)", S[20:70]),
        ("3' region (30-85)", S[30:]),
    ]
    
    centers = []
    labels = []
    
    for label, seq in sequences:
        if len(seq) >= window_size:
            result = analyze_promoter(seq, window_size)
            centers.append(result['center'])
            labels.append(label)
            print(f"  {label}: Center=({result['center'][0]:.2f}, {result['center'][1]:.2f})")
    
    plt.figure(figsize=(10, 6))
    centers_cg = [c[0] for c in centers]
    centers_ic = [c[1] for c in centers]
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(centers)))
    
    for i, (cg, ic) in enumerate(centers):
        plt.scatter([cg], [ic], s=400, marker='o', alpha=0.7, 
                   color=colors[i], edgecolors='black', linewidth=2, label=labels[i])
        plt.annotate(f'{i+1}', (cg, ic), ha='center', va='center', 
                    fontsize=12, fontweight='bold')
    
    plt.xlabel('C+G Content (%)', fontsize=12)
    plt.ylabel('Index of Coincidence (χ²)', fontsize=12)
    plt.title('Centers of Weight - Multiple Sequence Regions', fontsize=14, fontweight='bold')
    plt. legend(fontsize=9, loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('centers_pattern.png', dpi=300, bbox_inches='tight')
    
    print("Promoter Sequence Analysis")
    
    promoter = "GACACCATCGAATGGCGCAAAACCTTTCGCGGTATGGCATGATAGCGCCCGGAAGAGAGTCAATTCAGGGTGGTGAATGTGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTTTCCCGCGTGGTGAACCAG"
    
    print(f"  Promoter sequence length: {len(promoter)}")
    
    if len(promoter) >= window_size:
        promoter_results = analyze_promoter(promoter, window_size)
        p_center_cg, p_center_ic = promoter_results['center']
        
        print(f"  CG content (weighted): {promoter_results['avg_cg']:.2f}%")
        print(f"  IC (average): {promoter_results['avg_ic']:.2f}")
        print(f"  Center of weight: ({p_center_cg:.2f}, {p_center_ic:.2f})")
        print(f"  Number of windows: {promoter_results['window_count']}")
        
        fig3 = plot_pattern(promoter_results['cg_values'], promoter_results['ic_values'],
                            title=f"DNA Pattern - Promoter Sequence ({promoter_results['window_count']} windows)")
        plt.savefig('promoter_pattern.png', dpi=300, bbox_inches='tight')
    
    plt.show()
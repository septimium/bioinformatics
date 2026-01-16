# Use the json file to create an engine that is able to use the transition probabilities in order to produce an output similar to the orignal text from which 
# the training was made. One has to create two separate apps, one producing DNA sequences and one producing english text.
import random
import json

class DNASequenceGenerator:
    def __init__(self, transition_probs=None):
        self.nucleotides = ['A', 'T', 'G', 'C']
        
        if transition_probs is None:
            self.transition_probs = {
                'A': {'A': 0.25, 'T': 0.30, 'G': 0.25, 'C': 0.20},
                'T': {'A': 0.30, 'T': 0.20, 'G': 0.20, 'C': 0.30},
                'G': {'A': 0.20, 'T': 0.25, 'G': 0.30, 'C': 0.25},
                'C': {'A': 0.25, 'T': 0.30, 'G': 0.25, 'C': 0.20}
            }
        else:
            self.transition_probs = transition_probs
    
    def load_from_json(self, json_file):
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        if 'transition_probabilities' in data:
            self.transition_probs = data['transition_probabilities']
    
    def train_from_sequence(self, sequence):
 
        counts = {n: {m: 0 for m in self.nucleotides} for n in self.nucleotides}
        
        for i in range(len(sequence) - 1):
            current = sequence[i].upper()
            next_base = sequence[i + 1].upper()
            
            if current in self.nucleotides and next_base in self.nucleotides:
                counts[current][next_base] += 1
        
        self.transition_probs = {}
        for current in self.nucleotides:
            total = sum(counts[current].values())
            if total > 0:
                self.transition_probs[current] = {
                    next_base: counts[current][next_base] / total
                    for next_base in self.nucleotides
                }
            else:
                self.transition_probs[current] = {n: 0.25 for n in self.nucleotides}
    
    def get_next_nucleotide(self, current):
        if current not in self.transition_probs:
            return random.choice(self.nucleotides)
        
        transitions = self.transition_probs[current]
        nucleotides = list(transitions.keys())
        probabilities = list(transitions.values())
        
        return random.choices(nucleotides, weights=probabilities, k=1)[0]
    
    def generate_sequence(self, length=100, start_nucleotide=None):
        if start_nucleotide is None:
            current = random.choice(self.nucleotides)
        else:
            current = start_nucleotide.upper()
        
        sequence = [current]
        
        for _ in range(length - 1):
            next_nuc = self.get_next_nucleotide(current)
            sequence.append(next_nuc)
            current = next_nuc
        
        return ''.join(sequence)
    
    def format_sequence(self, sequence, line_length=60):
        lines = []
        for i in range(0, len(sequence), line_length):
            lines.append(sequence[i:i + line_length])
        return '\n'.join(lines)
    
    def analyze_sequence(self, sequence):
        sequence = sequence.upper()
        total = len(sequence)
        
        counts = {n: sequence.count(n) for n in self.nucleotides}
        percentages = {n: (counts[n] / total) * 100 for n in self.nucleotides}
        
        gc_content = ((counts['G'] + counts['C']) / total) * 100
        
        return {
            'length': total,
            'counts': counts,
            'percentages': percentages,
            'gc_content': gc_content
        }


def main():
    print("=" * 70)
    print("DNA Sequence Generator using Markov Chain")
    print("=" * 70)
    print()
    
    generator = DNASequenceGenerator()
    
    print("Training on sample DNA sequence...")
    sample_dna = "ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTACGTACGTACGTAGCTAGCTAGCT"
    generator.train_from_sequence(sample_dna)
    print("Training complete!")
    print()
    
    print("Generated DNA Sequences:")
    print("-" * 70)
    
    for i in range(3):
        print(f"\nSequence {i + 1}:")
        sequence = generator.generate_sequence(length=120)
        formatted = generator.format_sequence(sequence, line_length=60)
        print(formatted)
        
        analysis = generator.analyze_sequence(sequence)
        print(f"\nAnalysis:")
        print(f"  Length: {analysis['length']} bp")
        print(f"  A: {analysis['counts']['A']} ({analysis['percentages']['A']:.1f}%)")
        print(f"  T: {analysis['counts']['T']} ({analysis['percentages']['T']:.1f}%)")
        print(f"  G: {analysis['counts']['G']} ({analysis['percentages']['G']:.1f}%)")
        print(f"  C: {analysis['counts']['C']} ({analysis['percentages']['C']:.1f}%)")
        print(f"  GC Content: {analysis['gc_content']:.1f}%")
        print("-" * 70)
    
    
    while True:
        user_input = input("\nEnter sequence length to generate (or 'quit' to exit): ")
        
        if user_input.lower() == 'quit':
            break
        
        try:
            length = int(user_input)
            sequence = generator.generate_sequence(length=length)
            print("\nGenerated DNA Sequence:")
            print(generator.format_sequence(sequence))
            
            # Ask if user wants analysis
            if input("\nShow analysis? (y/n): ").lower() == 'y':
                analysis = generator.analyze_sequence(sequence)
                print(f"\nAnalysis:")
                print(f"  Length: {analysis['length']} bp")
                print(f"  GC Content: {analysis['gc_content']:.1f}%")
                for nuc in ['A', 'T', 'G', 'C']:
                    print(f"  {nuc}: {analysis['counts'][nuc]} ({analysis['percentages'][nuc]:.1f}%)")
        
        except ValueError:
            print("Please enter a valid number.")


if __name__ == "__main__":
    main()
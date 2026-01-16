# Use an english text of a length of 300 letters and calculate the transition probbilities between thw words inside the text.
# Store that as JSON file. Represent each new word with one unique symbol.
import json
import re
from collections import defaultdict

def clean_text(text):
    text = text.lower()
    text = re.sub(r'[^\w\s]', '', text)
    text = ' '.join(text.split())
    return text

def extract_words(text):
    cleaned = clean_text(text)
    words = cleaned.split()
    return words

def assign_symbols(words):
    unique_words = []
    seen = set()
    
    for word in words:
        if word not in seen:
            unique_words.append(word)
            seen.add(word)
    
    word_to_symbol = {word: f"S{i}" for i, word in enumerate(unique_words)}
    symbol_to_word = {f"S{i}": word for i, word in enumerate(unique_words)}
    
    return word_to_symbol, symbol_to_word

def calculate_word_transitions(words, word_to_symbol):
    transition_counts = defaultdict(lambda: defaultdict(int))
    
    for i in range(len(words) - 1):
        current_word = words[i]
        next_word = words[i + 1]
        current_symbol = word_to_symbol[current_word]
        next_symbol = word_to_symbol[next_word]
        transition_counts[current_symbol][next_symbol] += 1
    
    transition_probs = {}
    all_symbols = set(word_to_symbol.values())
    
    for from_symbol in all_symbols:
        total = sum(transition_counts[from_symbol].values())
        if total > 0:
            transition_probs[from_symbol] = {
                to_symbol: transition_counts[from_symbol][to_symbol] / total
                for to_symbol in all_symbols
                if transition_counts[from_symbol][to_symbol] > 0
            }
        else:
            transition_probs[from_symbol] = {}
    
    return transition_probs, transition_counts

def print_results(words, word_to_symbol, symbol_to_word, transition_probs, transition_counts):
    print("\n" + "="*70)
    print("WORD TRANSITION ANALYSIS")
    print("="*70)
    
    print(f"\nTotal words in text: {len(words)}")
    print(f"Unique words: {len(word_to_symbol)}")
    
    print("\n" + "-"*70)
    print("WORD-SYMBOL MAPPING:")
    print("-"*70)
    for word, symbol in sorted(word_to_symbol.items(), key=lambda x: x[1]):
        print(f"  {symbol}: '{word}'")
    
    print("\n" + "-"*70)
    print("TRANSITION PROBABILITIES:")
    print("-"*70)
    
    for from_symbol in sorted(transition_probs.keys()):
        from_word = symbol_to_word[from_symbol]
        print(f"\nFrom {from_symbol} ('{from_word}'):")
        
        if transition_probs[from_symbol]:
            sorted_transitions = sorted(
                transition_probs[from_symbol].items(),
                key=lambda x: x[1],
                reverse=True
            )
            
            for to_symbol, prob in sorted_transitions:
                to_word = symbol_to_word[to_symbol]
                count = transition_counts[from_symbol][to_symbol]
                print(f"  → {to_symbol} ('{to_word}'): {prob:.4f} ({prob*100:.2f}%) [count: {count}]")
        else:
            print(f"  (No transitions)")

def save_to_json(data, filename="assign3.json"):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    print(f"\n" + "="*70)
    print(f"Results saved to '{filename}'")
    print("="*70)

def main():
    print("Word Transition Probability Calculator")
    print("="*70)
    text = """
    The quick brown fox jumps over the lazy dog. The dog was sleeping under 
    a big tree. The fox was very clever and fast. It ran through the forest 
    and found a river. The river was clear and cold. The fox drank some water 
    and rested. After a while, the fox continued its journey through the forest. 
    The sun was shining bright in the sky.
    """
    
    print(f"\nOriginal text (length: {len(text)} characters):")
    print("-"*70)
    print(text.strip())
    print("-"*70)
    
    words = extract_words(text)
    
    word_to_symbol, symbol_to_word = assign_symbols(words)
    
    transition_probs, transition_counts = calculate_word_transitions(words, word_to_symbol)
    
    print_results(words, word_to_symbol, symbol_to_word, transition_probs, transition_counts)
    
    json_data = {
        "original_text": text.strip(),
        "text_length": len(text),
        "word_count": len(words),
        "unique_word_count": len(word_to_symbol),
        "word_sequence": words,
        "word_to_symbol_mapping": word_to_symbol,
        "symbol_to_word_mapping": symbol_to_word,
        "transition_probabilities": transition_probs,
        "transition_counts": {
            from_sym: dict(to_dict) 
            for from_sym, to_dict in transition_counts.items()
        },
        "description": "Word transition probabilities from English text analysis"
    }
    
    save_to_json(json_data)
    
    print("\nADDITIONAL STATISTICS:")
    print("-"*70)
    
    word_freq = defaultdict(int)
    for word in words:
        word_freq[word] += 1
    
    print("\nMost frequent words:")
    for word, count in sorted(word_freq.items(), key=lambda x: x[1], reverse=True)[:5]:
        symbol = word_to_symbol[word]
        percentage = (count / len(words)) * 100
        print(f"  {symbol} ('{word}'): {count} times ({percentage:.2f}%)")
    
    singleton_words = [word for word, count in word_freq.items() if count == 1]
    print(f"\nWords appearing only once: {len(singleton_words)}")

if __name__ == "__main__":
    main()
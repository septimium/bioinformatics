import json
import random

class EnglishTextGenerator:
    def __init__(self, json_file):
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        self.transition_probs = data['transition_probabilities']
        self.symbol_to_word = data['symbol_to_word_mapping']
        self.word_to_symbol = data['word_to_symbol_mapping']
        
    def get_next_symbol(self, current_symbol):
        if current_symbol not in self.transition_probs:
            return None
        
        transitions = self.transition_probs[current_symbol]
        
        if not transitions:
            return None
        
        symbols = list(transitions.keys())
        probabilities = list(transitions.values())
        
        next_symbol = random.choices(symbols, weights=probabilities, k=1)[0]
        return next_symbol
    
    def generate_text(self, num_words=50, start_word='the'):
        if start_word.lower() not in self.word_to_symbol:
            start_word = random.choice(list(self.word_to_symbol.keys()))
        
        current_symbol = self.word_to_symbol[start_word.lower()]
        generated_words = [self.symbol_to_word[current_symbol]]
        
        for _ in range(num_words - 1):
            next_symbol = self.get_next_symbol(current_symbol)
            
            if next_symbol is None:
                next_symbol = self.word_to_symbol['the']
            
            word = self.symbol_to_word[next_symbol]
            generated_words.append(word)
            current_symbol = next_symbol
        
        text = ' '.join(generated_words)
        text = text[0].upper() + text[1:]
        
        words_list = text.split()
        formatted_text = []
        for i, word in enumerate(words_list):
            formatted_text.append(word)
            if i > 0 and (i + 1) % random.randint(8, 12) == 0 and i < len(words_list) - 1:
                formatted_text[-1] += '.'
                if i + 1 < len(words_list):
                    words_list[i + 1] = words_list[i + 1].capitalize()
        
        formatted_text[-1] += '.'
        return ' '.join(formatted_text)


def main():
    generator = EnglishTextGenerator('assign3.json')
    
    print("English Text Generator")
    print()
    
    for i in range(3):
        print(f"Sample {i + 1}:")
        text = generator.generate_text(num_words=60, start_word='the')
        print(text)
        print()
    
    while True:
        user_input = input("\nEnter number of words to generate (or 'quit' to exit): ")
        
        if user_input.lower() == 'quit':
            break
        
        try:
            num_words = int(user_input)
            start = input("Enter starting word (press Enter for 'the'): ").strip()
            if not start:
                start = 'the'
            
            text = generator.generate_text(num_words=num_words, start_word=start)
            print("\nGenerated text:")
            print(text)
        except ValueError:
            print("Please enter a valid number.")


if __name__ == "__main__":
    main()
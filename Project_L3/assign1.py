# The melting temperature (Tm) is the temperature at which one-half of a particular DNA duplex will dissociate and become a single strand of DNA. 
# Primer length and sequence are of critical importance in designing the parameters of a successful amplification. 
# The melting temperature of a nucleic acid duplex increases both with its length, and with increasing GC content. 
# A simple formula for calculation of the (Tm) is:  Tm = 4(G + C) + 2(A + T) °C
# The actual Tm is influenced by the concentration of Mg2+ , K+ , and cosolvents. An alternative formula is:
# Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
# Implement an application that calculates the melting temperature of a DNA sequence using one of these formulas or both.
# Input = a string of DNA
# Output = temperature in celsius

import math
def formula1(dna):
    dna = dna.upper()
    g = dna.count('G')
    c = dna.count('C')
    a = dna.count('A')
    t = dna.count('T')
    tm = 4 * (g + c) + 2 * (a + t)
    return tm

def formula2(dna, na_conc=0.05):
    dna = dna.upper()
    length = len(dna)
    gc_percent = (dna.count('G') + dna.count('C')) / length * 100
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / length)
    return tm

if __name__ == "__main__":
    dna = input("Enter DNA sequence: ").strip()
    tm1 = formula1(dna)
    tm2 = formula2(dna)
    print(f"\nFirst formula: {tm1:.2f} °C")
    print(f"Second formula: {tm2:.2f} °C")

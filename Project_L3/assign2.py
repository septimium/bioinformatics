# Design an application that uses the sliding window method in order to read the melting temperature over the sequence S.
# Use a sliding window of 8 positions and choose a fasta file as input.
import math
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as file:
        for line in file:
            if not line.startswith(">"): 
                sequence += line.strip().upper()
    return sequence

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


def sliding_window_tm(seq, window_size=8):
    tm_form1 = []
    tm_form2= []
    positions = []

    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        tm_form1.append(formula1(window))
        tm_form2.append(formula2(window))
        positions.append(i + 1)

    return positions, tm_form1, tm_form2


def smooth_data(values, window_length=5, polyorder=2):
    if len(values) < window_length:
        return values  #
    if window_length % 2 == 0:
        window_length += 1
    return savgol_filter(values, window_length, polyorder)

def main():
    fasta_file = input("Enter FASTA filename: ").strip()
    seq = read_fasta(fasta_file)
    window_size = 8
    positions, tm_form1, tm_form2 = sliding_window_tm(seq, window_size)
    data1 = smooth_data(tm_form1, window_length=15, polyorder=2)
    data2 = smooth_data(tm_form2, window_length=15, polyorder=2)
    plt.figure(figsize=(10, 6))
    plt.plot(positions, data1, label="Formula 1", color="blue", linewidth=2)
    plt.plot(positions, data2, label="Formula 2", color="red", linewidth=2)
    plt.title(f"Melting Temperature (Tm) Profile\nWindow size = {window_size}")
    plt.xlabel("Start position of window")
    plt.ylabel("Tm (°C)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
# 1. Make an artificial DNA sequence of 200-400b in length, in which to simulate 3-4 transposable elements.
# 2. Implement a software application to detect the positions of these transposable elements (start, end) 
# within the created DNA sequence.

import random

def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

genome_length = random.randint(200, 400)
genome = list(random_dna(genome_length))

num_tes = random.randint(3, 4)

transposons = [random_dna(random.randint(15, 25)) for _ in range(num_tes)]

insert_positions = []

pos = random.randint(0, len(genome) - 1)
insert_positions.append(pos)
genome[pos:pos] = list(transposons[0])

pos2 = pos + random.randint(1, len(transposons[0]) - 5)
insert_positions.append(pos2)
genome[pos2:pos2] = list(transposons[1])

if num_tes >= 3:
    pos3 = pos + random.randint(3, len(transposons[0]) - 8)
    insert_positions.append(pos3)
    genome[pos3:pos3] = list(transposons[2])

if num_tes == 4:
    pos4 = random.randint(0, len(genome) - 1)
    insert_positions.append(pos4)
    genome[pos4:pos4] = list(transposons[3])

genome = "".join(genome)

detections = []

for te in transposons:
    start = 0
    while True:
        idx = genome.find(te, start)
        if idx == -1:
            break
        detections.append((te, idx, idx + len(te)))
        start = idx + 1  


print("\n=== GENOME LENGTH AFTER INSERTION ===")
print(len(genome))

print("\n=== TRANSPOSONS (SOME OVERLAPPING) ===")
for i, (te, pos) in enumerate(zip(transposons, insert_positions), 1):
    print(f"TE{i}: {te}  | Inserted at index: {pos}")

print("\n=== DETECTED TRANSPOSON POSITIONS ===")
for te, start, end in detections:
    print(f"Detected TE '{te}' from {start} to {end}")
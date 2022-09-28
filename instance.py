from math import ceil, inf
from collections import defaultdict
import numpy as np

# init constants
NUCLEOTIDES = np.array(['G', 'T', 'C', 'A'])

# get parameters
DNA_SIZE = int(input("Enter DNA size: "))
OLIGO_SIZE = int(input("Enter oligonucleotide size: "))
ERROR_PERCENTAGE = int(input("Enter error %: "))

# generate DNA sequence
DNA_SEQUENCE = "".join(np.random.choice(NUCLEOTIDES, size=DNA_SIZE))
FIRST_OLIGO = DNA_SEQUENCE[0:OLIGO_SIZE]

# generate ideal DNA spectrum
ideal_spectrum = defaultdict(int)
SPECTRUM_SIZE = DNA_SIZE - OLIGO_SIZE + 1

for i in range(SPECTRUM_SIZE):
    oligo = DNA_SEQUENCE[i:i+OLIGO_SIZE]
    ideal_spectrum[oligo] += 1

oligos = list(ideal_spectrum.keys())

NUM_ERRORS = ceil(SPECTRUM_SIZE * ERROR_PERCENTAGE / 100)

# simulate negative errors
removed_oligos = []

for _ in range(NUM_ERRORS):
    index = np.random.randint(1, len(oligos))
    removed = oligos.pop(index)
    ideal_spectrum.pop(removed)
    removed_oligos.append(removed)

# simulate positive errors
added_oligos = []
for _ in range(NUM_ERRORS):
    while (new_oligo := "".join(np.random.choice(NUCLEOTIDES, size=OLIGO_SIZE))) in oligos:
        continue

    added_oligos.append(new_oligo)
    ideal_spectrum[new_oligo] += 1

# write to file
file = open("data.txt", "w", encoding="utf8")
file.write(DNA_SEQUENCE + '\n')

file.write(F"{FIRST_OLIGO} {ideal_spectrum[FIRST_OLIGO]}\n")

ideal_spectrum.pop(FIRST_OLIGO)
final_oligos = list(ideal_spectrum.keys())

while len(final_oligos) > 0:
    index = np.random.randint(len(final_oligos))
    oligo = final_oligos.pop(index)
    rep_count = ideal_spectrum[oligo]
    rep_count = inf if rep_count > 2 else rep_count
    file.write(f"{oligo} {rep_count}\n")

file.close()

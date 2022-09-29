# pylint: disable=invalid-name, missing-module-docstring, missing-function-docstring, missing-class-docstring
from dataclasses import dataclass
import numpy as np
from math import inf

@dataclass
class Instance:
    original_sequence: str
    oligos: list[tuple[str, int or np.inf]]


def load_instance(file_path):
    with open(file_path, encoding='utf-8') as f:
        oligos = []
        dna_sequence = f.readline().strip()
        for line in f:
            split_line = line.split()
            oligos.append(
                (
                    split_line[0],
                    float(split_line[1]) if split_line[1] == "inf" else int(split_line[1])
                )
            )
        return Instance(dna_sequence, oligos)



def select_one(population):
    maximum = sum(list(population))
    selection_probs = [c/maximum for c in population]
    return population[np.random.choice(len(population), p=selection_probs)]


def calculateDistance(oligo1: str, oligo2: str):
    rowSize = len(oligo1)
    columnSize = len(oligo2)
    if rowSize != columnSize:
        return None

    for i in range(0, rowSize):
        if oligo1[i:rowSize] == oligo2[0: rowSize - i]:
            return i
    return 0


def assembleDNA(solution: list, oligo_size: int):
    result = solution[0]
    while len(solution) != 1:
        dist = calculateDistance(solution[0], solution[1])
        tmp = solution[1][oligo_size - dist: oligo_size]
        result += tmp
        solution.remove(solution[0])
    return result


def levenshteinDistance(string1 : str, string2 : str):
    """
    Implementation of http://web.archive.org/web/20120526085419/http://www.merriampark.com/ldjava.htm
    """

    n = len(string1)
    m = len(string2)

    if  m == 0:
        return n
    elif n == 0:
        return m

    prevCost = []
    cost = []

    for i in range(n+1):
        prevCost.append(i)
        cost.append(0)

    for j in range(1, m + 1):
        string2_char = string2[j - 1]
        cost[0] = j
        for i in range(1, n+1):
            value = 0 if string1[i-1] == string2_char else 1
            cost[i] = min(cost[i - 1] + 1, prevCost[i] + 1, prevCost[i - 1] + value)
        prevCost, cost = cost, prevCost

    return prevCost[n]

# def levenshteinDistance(token1, token2):
#     distances = np.zeros((len(token1) + 1, len(token2) + 1))

#     for t1 in range(len(token1) + 1):
#         distances[t1][0] = t1

#     for t2 in range(len(token2) + 1):
#         distances[0][t2] = t2
        
#     a = 0
#     b = 0
#     c = 0
    
#     for t1 in range(1, len(token1) + 1):
#         for t2 in range(1, len(token2) + 1):
#             if (token1[t1-1] == token2[t2-1]):
#                 distances[t1][t2] = distances[t1 - 1][t2 - 1]
#             else:
#                 a = distances[t1][t2 - 1]
#                 b = distances[t1 - 1][t2]
#                 c = distances[t1 - 1][t2 - 1]
                
#                 if (a <= b and a <= c):
#                     distances[t1][t2] = a + 1
#                 elif (b <= a and b <= c):
#                     distances[t1][t2] = b + 1
#                 else:
#                     distances[t1][t2] = c + 1

#     return distances[len(token1)][len(token2)]


if __name__ == "__main__":
    print(
        levenshteinDistance(
            "GTTGCAAATATTCTTGTCGGGGAACGCTCTTTAGCGTCTCTCCATGTAGGAGGAGAGCACACACCCTCCTAAGGCATGTCTTACTCCCATATATGTACAC",
            'GTTGCAAATATTCTTGTCGGGGAACGCTCTTTAGCGTCTCTCCATGTAGGAGGAGAGCACACACCCTCCTAAGGCATGTCTTACTCCCATATATGTACAC'
        )
    )
    # print(calculateDistance('ATGT', 'ATGT'))
    # print(calculateDistance('ATGT', 'TGTA'))
    # print(calculateDistance('ATGT', 'GTAT'))
    # print(calculateDistance('ATGT', 'TGAT'))
    # print(calculateDistance('ATGT', 'GAAA'))
    # print(calculateDistance('TTATGAT', 'TATGATG'))
    # print(assembleDNA(['ATGC', 'GCTC', 'TCTA', 'ATGC'], 4))
    # print("module test suite")
    # instance = load_instance("data.txt")
    # print(instance.original_sequence)
    # print(instance.oligos)

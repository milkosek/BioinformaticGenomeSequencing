from collections import defaultdict
import numpy as np


def select_one(population):
    maximum = sum(list(population))
    selection_probs = [c/maximum for c in population]
    return population[np.random.choice(len(population), p=selection_probs)]


if __name__ == "__main__":
    print("module test suite")
    population = [0.1, 0.2, 0.3, 0.4]
    results = defaultdict(int)
    for _ in range(100000):
        results[select_one(population)] += 1
    print(results)

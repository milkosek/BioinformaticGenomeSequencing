from collections import defaultdict
import numpy as np
from utils import Instance, assembleDNA, calculateDistance, levenshteinDistance, load_instance


class AntColony:
    def __init__(
        self,
        inst: Instance,
        ants_count: int = 50,
        alpha: float = 2,
        beta: float = 3,
        evaporation: float = 0.2,
        iterations: int = 100
        ):
        self.instance = inst
        self.ants_count = ants_count
        self.alpha = alpha
        self.beta = beta
        self.evaporation = evaporation
        self.iterations = iterations
        self.best_levenshtein = float('inf')
        self.best_solution = ""
        self.best_route = []

        self.oligos_count = len(instance.oligos)
        self.oligo_size = len(instance.oligos[0][0])
        self.dna_size = len(instance.original_sequence)

        self.oligos_map = defaultdict(int)
        self.oligos_max_count = defaultdict(int)
        self.oligos_optimal_successors = defaultdict(list)
        self.oligos_list = [oligo[0] for oligo in self.instance.oligos]
        self.pheromone_matrix = None
        self.distance_matrix = None
        self._init_structures()


    def _map_oligos_to_indices(self):
        count = 0
        for oligo in self.instance.oligos:
            self.oligos_map[oligo[0]] = count
            self.oligos_max_count[oligo[0]] = oligo[1]
            count += 1

    def _get_optimal_successors(self):
        for oligo in self.oligos_list:
            current = self.oligos_map[oligo]
            self.oligos_optimal_successors[oligo] = list(
                filter(
                    lambda successor: 0 < self.distance_matrix[current][self.oligos_map[successor]] < 4,
                    self.oligos_list
                )
            )


    def _generate_pheromone_matrix(self):
        self.pheromone_matrix = [
            [0.5 for _ in range(self.oligos_count)]
            for _ in range(self.oligos_count)
        ]


    def _generate_distance_matrix(self):
        size = self.oligos_count
        self.distance_matrix = [
            [0 for _ in range(size)]
            for _ in range(size)
        ]
        for i in range(size):
            for j in range(size):
                self.distance_matrix[i][j] = calculateDistance(
                        self.instance.oligos[i][0],
                        self.instance.oligos[j][0]
                    )


    def _init_structures(self):
        self._map_oligos_to_indices()
        self._generate_distance_matrix()
        self._generate_pheromone_matrix()
        self._get_optimal_successors()


    # def _update_pheromone(self, pheromones_matrices: list):
    #     # evaporation
    #     for i in range(self.oligos_count):
    #         for j in range(self.oligos_count):
    #             self.pheromone_matrix[i][j] *= (1 - self.evaporation)
    #     # pheromone from ants
    #     for pheromone_matrix in pheromones_matrices:
    #         for i in range(self.oligos_count):
    #             for j in range(self.oligos_count):
    #                 self.pheromone_matrix[i][j] += pheromone_matrix[i][j]
    #                 if self.pheromone_matrix[i][j] > 100:
    #                     self.pheromone_matrix[i][j] = 100

    def _update_pheromone(self):
        # evaporation
        for i in range(self.oligos_count):
            for j in range(self.oligos_count):
                self.pheromone_matrix[i][j] *= (1 - self.evaporation)
        # pheromone from current best solution
        for i in range(len(self.best_route) - 1):
            idx1 = self.oligos_map[self.oligos_list[i]]
            idx2 = self.oligos_map[self.oligos_list[i + 1]]
            self.pheromone_matrix[idx1][idx2] += (20 / self.best_levenshtein)
            


    def _update_best_solution(self, new_route, new_dna, new_levenshtein):
        if new_levenshtein < self.best_levenshtein:
            self.best_route = new_route
            self.best_solution = new_dna
            self.best_levenshtein = new_levenshtein



    def _ant_get_weight(self, node1, node2):
        idx1 = self.oligos_map[node1]
        idx2 = self.oligos_map[node2]
        pheromone = self.pheromone_matrix[idx1][idx2]
        distance = self.distance_matrix[idx1][idx2]
        weight = (pheromone**self.alpha) * ((1 / distance)**self.beta)
        return weight


    def _ant_run(self):
        # initialize local pheromone accumulator
        # ant_pheromone_matrix = [
        #     [float(0) for _ in range(self.oligos_count)]
        #     for _ in range(self.oligos_count)
        # ]
        # local counter of nodes used in solution
        nodes_use_count = defaultdict(int)
        for oligo in self.oligos_list:
            nodes_use_count[oligo] = 0

        route = [self.oligos_list[0]]
        solution_size = self.oligo_size
        current_node = route[0]
        nodes_use_count[current_node] += 1
        used_nodes = []
        if nodes_use_count[current_node] == self.oligos_max_count[current_node]:
            used_nodes.append(current_node)

        while solution_size < self.dna_size:
            # get next possible oligos
            # (filter oligos with max number of uses and solution overflow)
            successors = list(filter(
                lambda oligo: oligo not in used_nodes and
                    self.distance_matrix[self.oligos_map[current_node]][self.oligos_map[oligo]] + solution_size <= self.dna_size
                ,
                self.oligos_optimal_successors[current_node]))
            # print(f"successors: {successors}")
            # there is no successor, finish solution build
            if len(successors) == 0:
                break
            # calculate weights for random selection
            weights = [
                self._ant_get_weight(current_node, successor)
                for successor in successors
            ]
            # max_w = 0
            # min_w = float('inf')
            # for value in weights:
            #         if value > max_w:
            #             max_w = value
            #         if value < min_w:
            #             min_w = value
            # print(f"weights: [{min_w}, {max_w}]")
            # get successor
            [next_node] = np.random.choice(successors, 1, weights)
            # update local solution
            route.append(next_node)
            nodes_use_count[next_node] += 1
            solution_size += self.distance_matrix[
                self.oligos_map[current_node]
            ][
                self.oligos_map[next_node]
            ]
            if nodes_use_count[next_node] == self.oligos_max_count[next_node]:
                used_nodes.append(current_node)
            current_node = next_node
        
        # distribute pheromone based on route
        generated_dna = assembleDNA(list(route), self.oligo_size)
        route_rank = levenshteinDistance(self.instance.original_sequence, generated_dna)
        self._update_best_solution(route, generated_dna, route_rank)
        # print(route_rank)
        # skip next calculations if best solution was found
        if route_rank == 0:
            return route, generated_dna, True
            # return route, generated_dna, ant_pheromone_matrix, True
        # for i in range(len(route) - 1):
        #     idx1 = self.oligos_map[route[i]]
        #     idx2 = self.oligos_map[route[i + 1]]
        #     ant_pheromone_matrix[idx1][idx2] += (1.0 / route_rank)

        # print(ant_pheromone_matrix)
        # return route, generated_dna, ant_pheromone_matrix, False
        return route, generated_dna, False


    def run(self):
        for iteration in range(self.iterations):
            stop = False
            print(f"Starting iteration #{iteration}")
            # pheromones_matrices = []
            for _ in range(self.ants_count):
                # print(f"Ant #{ant} generates solution")
                ant_route, ant_dna, found_best = self._ant_run()
                # ant_route, ant_dna, ant_pheromone_matrix, found_best = self._ant_run()
                if found_best:
                    stop = True
                    break
                # pheromones_matrices.append(ant_pheromone_matrix)

            if stop:
                break
            self._update_pheromone()
            print(self.best_levenshtein)
            # max_pheromone = 0
            # min_pheromone = float('inf')
            # for row in self.pheromone_matrix:
            #     for value in row:
            #         if value > max_pheromone:
            #             max_pheromone = value
            #         if value < min_pheromone:
            #             min_pheromone = value
            # print(f"min: {min_pheromone} max: {max_pheromone}")
            # self._update_pheromone(pheromones_matrices)
        return self.best_solution



if __name__ == "__main__":
    instance = load_instance("data.txt")
    alpha = 2
    beta = 3
    ants_count = 20
    evaporation = 0.1
    iterations = 50
    ac_instance = AntColony(
        instance,
        ants_count=ants_count,
        alpha=alpha,
        beta=beta,
        evaporation=evaporation,
        iterations=iterations
    )
    # for oligo in ac_instance.oligos_optimal_successors:
    #     print(oligo)
    #     print(ac_instance.oligos_optimal_successors[oligo])
    instance_solution = ac_instance.run()
    print(f"Original sequence:\n{instance.original_sequence} {len(instance.original_sequence)}")
    print(f"Generated sequence:\n{instance_solution} {len(instance_solution)}")
    print(f"Levenshtein distance: {levenshteinDistance(instance.original_sequence, instance_solution)}")
    
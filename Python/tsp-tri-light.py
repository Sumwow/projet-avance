"""
Algorithme génétique pour le TSP
Vincent Dugat 2025
"""
import random,math
import tsplib95 # fonctions pour lire les fichiers TSPLIB et faire des calculs
import pyproj # pour les coordonnées géométriques
import matplotlib.pyplot as plt # pour les graphiques
from tsp_tools import *
from tsp_ga_fct import *

# variables globales : paramètres de l'algorithme génétique

population_size = 30
generations = 1000
mutation_rate = 0.10
tournament_size = int(0.5*population_size)

file = "att48"
instance = load_tsp_file(file+".tsp")
#instance = load_tsp_file("dja1436.tsp")

cities = list(instance.get_nodes())
instance_dict = instance.as_name_dict()
#print(instance)
coord = instance_dict['node_coords']
edge_type = instance_dict['edge_weight_type']

distance_fct,graphique = select_fct(edge_type)

# début des calculs
def main():

    population = random_population(population_size, cities)
    best_individual = population[0] #
    for generation in range(generations):
        selected = tournament_selection(population, fitness, tournament_size,distance_fct, coord)
        offspring = []
        for i in range(0, population_size, 2):
            child_a = ordered_crossover(selected[i], selected[i + 1])
            child_b = ordered_crossover(selected[i + 1], selected[i])
            offspring.append(child_a)
            offspring.append(child_b)

        for child in offspring:
            swap_mutation(child, mutation_rate)
        offspring = sorted(offspring, key=lambda indiv: fitness(indiv,distance_fct, coord)) # on trie
        population = offspring
        worste = max(population,key=lambda indiv: fitness(indiv,distance_fct,coord)) # la pire
        ind = population.index(worste)
        population[ind] = random.sample(cities, len(cities)) # on la remplace par du neuf
        best_generation = min(population, key=lambda indiv:fitness(indiv,distance_fct, coord)) #la meilleure
        worste_generation = max(population,key=lambda indiv: fitness(indiv,distance_fct,coord)) # la pire le retour
        ind = population.index(worste_generation)
        if fitness(best_generation,distance_fct, coord) < fitness(best_individual,distance_fct, coord): # mémoire générationnelle
            best_individual = best_generation
        population[ind] = best_individual # on le remet dans la piscine

    print_solution(best_individual,fitness,distance_fct,coord,file,instance,graphique)

main()

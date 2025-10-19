import random,math
import subprocess
import time
import os
import sys
import tsplib95 # fonctions pour lire les fichiers TSPLIB et faire des calculs
import pyproj # pour les coordonnées géométriques
import matplotlib.pyplot as plt # pour les graphiques
from tsp_tools import *

def random_population(n, cities):
    # crée une tournée aléatoire
    population = []
    for _ in range(n):
        individual = random.sample(cities, len(cities))
        population.append(individual)
    #print(population)
    return population

def tournament_selection(population, fitness_function, tournament_size,distance_fct,coord):
    selected = []
    for _ in range(len(population)):
        competitors = random.sample(population, tournament_size)# on sélectionne tournament_size individus
        winner = min(competitors,key=lambda indiv:fitness(indiv,distance_fct, coord)) # on prend le meilleur
        selected.append(winner)
    return selected # nouvelle population

def ordered_crossover(parent_a, parent_b):
    # croisement génétique
    start, end = sorted(random.sample(range(len(parent_a)), 2))
    child = [None] * len(parent_a)
    child[start:end] = parent_a[start:end]

    for gene in parent_b:
        if gene not in child:
            i = child.index(None)
            child[i] = gene

    return child

def swap_mutation(individual, mutation_rate):
    # mutation au hasard
    for i in range(len(individual)):
        if random.random() < mutation_rate:
            j = random.randint(0, len(individual) - 1)
            individual[i], individual[j] = individual[j], individual[i]

def print_solution(best_individual,fitness,distance_fct,coord,file,instance,graphique):
    dist = fitness(best_individual,distance_fct, coord)
    print(f"Best GA solution: {best_individual} longueur = {dist}")
    graphique(best_individual,"calcul",coord)
    #print(f"Longueur du tour canonique : {instance.trace_canonical_tour()}")
    if os.path.isfile(file+'.opt.tour'):
            opt = tsplib95.load(file+'.opt.tour')
            print(f"Tour optimum : {opt.tours[0]}, longueur : {fitness(opt.tours[0],distance_fct,coord)}")
            graphique(opt.tours[0],"optimal",coord)
    else:
        print("Pas de fichier solution présent")

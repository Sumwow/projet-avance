import random,math
import subprocess
import time
import os
import sys
import tsplib95 # fonctions pour lire les fichiers TSPLIB et faire des calculs
import pyproj # pour les coordonnées géométriques
import matplotlib.pyplot as plt # pour les graphiques

# constantes pour les coordonnées géographiques
RRR = 6378.388
PI = 3.141592

def load_tsp_file(filename):
    print("Lecture fichier :",filename)
    problem = tsplib95.loaders.load(filename)
    return problem

def show_tour(tour_list,com,coord):
    x = [coord[i][0] for i in tour_list] # attention à la numérotation dans tour_list
    y = [coord[i][1] for i in tour_list]
    x.append(x[0])
    y.append(y[0])
    plt.plot(x,y)
    plt.xlabel(com)
    plt.show()
    plt.close()

def dist_mat(i,j):
    if i<=j:
        i,j = j,i
    return dist_matrix[i-1][j-1]

def distance_eucl(city_a, city_b,coord)->int:
    # Calculate the distance between city_a and city_b
    # distance euclidienne
    x_a = coord[city_a][0]
    y_a = coord[city_a][1]
    x_b = coord[city_b][0]
    y_b = coord[city_b][1]
    x_offset = x_b - x_a
    y_offset = y_b - y_a
    dist = math.sqrt(x_offset*x_offset+y_offset*y_offset)
    return dist

def distance_att(city_a, city_b,coord):
    x_a = coord[city_a][0]
    y_a = coord[city_a][1]
    x_b = coord[city_b][0]
    y_b = coord[city_b][1]
    x_offset = x_b - x_a
    y_offset = y_b - y_a
    dist = math.sqrt((x_offset*x_offset+y_offset*y_offset)/10.0)
    tij = (int(dist+0.5));
    if (tij<dist):
        dist = tij + 1;
    else:
        dist = tij;
    #print(f">>>> {city_a} {city_b} {x_a} {y_a} {x_b} {y_b} {dist}")
    return dist

P = pyproj.Proj(proj='utm', zone=31, ellps='WGS84', preserve_units=True)
G = pyproj.Geod(ellps='WGS84')

def LonLat_To_XY(Lon, Lat):
    return P(Lon, Lat)

def XY_To_LonLat(x,y):
    return P(x, y, inverse=True)

def distance(Lat1, Lon1, Lat2, Lon2):
    return G.inv(Lon1, Lat1, Lon2, Lat2)[2]

def latitude(city,coord):
        deg = int( coord[city][0] )
        min = coord[city][0] - deg
        return PI * (deg + 5.0 * min / 3.0 ) / 180.0

def longitude(city,coord):
        deg = int( coord[city][1] )
        min = coord[city][1] - deg
        return PI * (deg + 5.0 * min / 3.0 ) / 180.0

def dist_geo2km(city_a,city_b,coord):
    # (lat,long)
    q1 = math.cos( longitude(city_a,coord) - longitude(city_b,coord) )
    q2 = math.cos( latitude(city_a,coord) - latitude(city_b,coord) )
    q3 = math.cos( latitude(city_a,coord) + latitude(city_b,coord) )
    dab = (int) ( RRR * math.acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) + 1.0)
    return dab

def fitness(individual,distance_fct,coord):
# calcul de la longueur de la tournée avec la distance passée en paramètre
    total_distance = 0
    for i in range(len(individual) - 1):
        total_distance += distance_fct(individual[i], individual[i + 1],coord)
    total_distance += distance_fct(individual[-1], individual[0],coord)  # Return to the first city
    return total_distance

def show_tour_geo(tour_list,com,coord):
    x,y = [],[]
    for city in tour_list:
        lg,lat = longitude(city,coord),latitude(city,coord)
        X,Y = LonLat_To_XY(lg,lat)
        x.append(X)
        y.append(Y)
    x.append(x[0])
    y.append(y[0])
    plt.plot(x,y)
    plt.xlabel(com)
    plt.show()
    plt.close()

def valid(tour):
    # vérifie que chaque sommet n'est présent qu'une fois
    return len([num for num in tour if tour.count(num) > 1])

def select_fct(edge_type):
    match edge_type: #choix de la fonction distance pour fitness
        case 'GEO':
            distance_fct = dist_geo2km
            graphique = show_tour_geo
            print("distance géométrique")
        case 'EXPLICIT':
            distance_fct = dist_mat
            #graphique = pas prévu
        case 'EUC_2D':
            distance_fct = distance_eucl
            graphique = show_tour
            print("distance euclidienne")
        case 'ATT':
            distance_fct = distance_att
            graphique = show_tour
            print("distance ATT")
        case _:
            print('Inconnu ou non géré')
            exit()
    return distance_fct,graphique

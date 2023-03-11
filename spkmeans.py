import numpy as np
import pandas as pd
import math
import sys
import spkmeansmodule as cmod

"""
Func implements the k-means++ centroid initialization algorythem.
@param k: The ammount of centroids to be chosen.
@param points: A matrix containing all of the points and their indices.
"""
def choose_centroids(k: int, points: np.ndarray) -> np.ndarray: #will be changed
    centroids = np.ndarray(shape = (k, len(points[0])))
    selected = np.random.choice(len(points))
    centroids[0] = points[selected]
    points = np.delete(points, selected, axis = 0)
    min_dists = np.ndarray(shape = len(points))
    for i in range(1, k):
        if i < k:
            update_min_dists(points, centroids[i - 1], min_dists, i == 1)
        selected = np.random.choice(np.asarray([i for i in range(len(points))]), p = calc_prob(min_dists))
        centroids[i] = points[selected]
        points = np.delete(points, selected, axis = 0)
        min_dists = np.delete(min_dists, selected)
    return centroids

"""
Func updates the minimum distance from the centroids after intializing a new centroid.
@param points: The remaining pointswhoch the new centroid will be chosen from.
@param new_centroid: The last selected centroid.
@param min_dists: An array containgin the current minimal distances od each point from the centroids.
@param first: Indicates wether min_dists has been initialized.
"""
def update_min_dists(points: np.ndarray, new_centroid: np.ndarray, min_dists: np.ndarray, first: bool) -> None:
    for i in range(len(points)):
        if first:
             min_dists[i] = distance(points[i], new_centroid)
        else:    
            min_dists[i] = min(min_dists[i], distance(points[i], new_centroid))

"""
Func computes the distance between 2 data points, ignoring the first (index) coloumn.
@param point1: The first point.
@param point2: The second point.
@returns: The distance between the points.
"""
def distance(point1: np.ndarray, point2: np.ndarray) -> float:
    return math.sqrt(sum([math.pow(point1[i] - point2[i], 2) for i in range(1, point1.size)]))

"""
Func computes the weighted probability of each one of the elements.
@param min_dists: An array containing the minimal distances of each point from the centroids.
@returns: The probabilites - min_dist/sum_min_dist
"""
def calc_prob(min_dists: np.ndarray) -> np.ndarray:
    probs = min_dists.copy()
    sum_min_dists = min_dists.sum()
    for i in range(len(probs)):
        probs[i] /= sum_min_dists
    return probs

def print_mat(points: list) -> None:
    for point in points:
        print_line(point)

def print_line(point: list) -> None:
    line = ""
    for cor in point:
        if type(cor) == int:
            line += '%d' % cor
        else:
            line += '%.4f' % cor
        line += ','
    print(line[:-1])
    
def main(k: int, goal: str, path: str, hur: bool): #will be changed
    points = pd.read_csv(path, header = None).to_numpy()
    if goal == "spk":
        np.random.seed(0)
        output = cmod.spk()
    elif goal == "wam":
        output = cmod.wam(points)
    elif goal == "ddg":
        output = cmod.ddg(points)
    elif goal == "gl":
        output = cmod.gl(points)
    else:
        output = cmod.jacobi(points)


if __name__ == "__main__":
    goals = ["spk", "wam", "ddg", "gl", "jacobi"]
    argv_len = len(sys.argv)
    if argv_len < 3 or argv_len > 4:
        print("Invalid input")
    elif argv_len == 4:
        if sys.argv[1].isdecimal():
            if sys.argv[2] not in goals:
                print("Invalid goal!")
            else:
                main(int(sys.argv[1]), sys.argv[2], sys.argv[3], True)
        else:
            print("Invalid number of clusters")
    else:
        if sys.argv[1] not in goals:
            print("Invalid goal!")
        else:
            main(0, sys.argv[1], sys.argv[2], False)

import numpy as np



def tanimoto(x, y):
    c = np.sum(x*y)
    a = np.sum(x)
    b = np.sum(y)
    
    return c/(a+b-c)


def euclidean_dist(x, y):
    return np.linalg.norm(x-y)


def cosine_sim(x, y):
    return np.dot(x, y)/(np.linalg.norm(x)*np.linalg.norm(y))
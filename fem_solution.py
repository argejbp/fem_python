import numpy as np 
import scipy as sc

def steady_state_solution(KG, FG):
    sol = np.linalg.inv(KG).dot(FG)
    return sol



import numpy as np 
import scipy as sc

def steady_state_solution(KG, FG, active_nodes):
    sol = np.linealg.solve(KG[active_nodes[:, None], active_nodes], FG[active_nodes, 0])
    return sol

def modal_solver(KG, MG, active_nodes):
    MGi = np.linealg.inv(MG[active_nodes[:, None], active_nodes])
    A = KG[active_nodes[:, None], active_nodes]*MGi
    v, d = sc.linalg.eig(A)
    return v, d

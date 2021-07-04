import numpy as np
from numpy.linalg import multi_dot
import math

pi = 4*math.atan(1)

#Heat Transfer Elements

#Linear Triangle Element
def Triangle_Element_HT(k, e, IN, Coordinates, t):
    pos_i = int(IN[0]) - 1
    pos_j = int(IN[1]) - 1
    pos_k = int(IN[2]) - 1

    xi = Coordinates[pos_i, 0]
    yi = Coordinates[pos_i, 1]
    xj = Coordinates[pos_j, 0]
    yj = Coordinates[pos_j, 1]
    xk = Coordinates[pos_k, 0]
    yk = Coordinates[pos_k, 1]


    A = 1/2*((xi*yi - xj*yi) + (xk*yi - xi*yk) + (xj*yk - xk*yj))
    
    beta1 = yj - yk
    beta2 = yk - yi
    beta3 = yi - yj

    gamma1 = xk - xj 
    gamma2 = xi - xk
    gamma3 = xj - xi

    Ke = k*t/(4*A)*np.array([[beta1**2 + gamma1**2, beta1*beta2 + gamma1*gamma2, beta1*beta3 + gamma1*gamma3], [beta1*beta2 + gamma1*gamma2, beta2**2 + gamma2**2, beta2*beta3 + gamma2*gamma3], [beta1*beta3 + gamma1*gamma3, beta2*beta3 + gamma2*gamma3, beta3**2 + gamma3**2]])
    return Ke

def Rectangular_Element_HT(prop, e, IN, Coordinates):
	#Coming soon :)
    pass

#Linear Triangle Element (axisymmetric element)
def Axisym_Triangle_Element_HT(k, e, IN, Coordinates):
    pos_i = int(IN[0]) - 1
    pos_j = int(IN[1]) - 1
    pos_k = int(IN[2]) - 1

    ri = Coordinates[pos_i, 0]
    zi = Coordinates[pos_i, 1]
    rj = Coordinates[pos_j, 0]
    zj = Coordinates[pos_j, 1]
    rk = Coordinates[pos_k, 0]
    zk = Coordinates[pos_k, 1]

    r_ = (ri + rj + rk)/3

    A = 1/2*((ri*zi - rj*zi) + (rk*zi - ri*zk) + (rj*zk - rk*zj))
    
    beta1 = zj - zk
    beta2 = zk - zi
    beta3 = zi - zj

    gamma1 = rk - rj 
    gamma2 = ri - rk
    gamma3 = rj - ri

    Ke = pi*k*r_/(2*A)*np.array([[beta1**2 + gamma1**2, beta1*beta2 + gamma1*gamma2, beta1*beta3 + gamma1*gamma3], [beta1*beta2 + gamma1*gamma2, beta2**2 + gamma2**2, beta2*beta3 + gamma2*gamma3], [beta1*beta3 + gamma1*gamma3, beta2*beta3 + gamma2*gamma3, beta3**2 + gamma3**2]])
    return Ke

#Solid Mechanics Elements

#Linear Triangle Element
def Triangle_Element_SM(E, v, e, IN, Coordinates, t):
    pos_i = int(IN[0]) - 1
    pos_j = int(IN[1]) - 1
    pos_k = int(IN[2]) - 1

    xi = Coordinates[pos_i, 0]
    yi = Coordinates[pos_i, 1]
    xj = Coordinates[pos_j, 0]
    yj = Coordinates[pos_j, 1]
    xk = Coordinates[pos_k, 0]
    yk = Coordinates[pos_k, 1]

    A = 1/2*((xi*yi - xj*yi) + (xk*yi - xi*yk) + (xj*yk - xk*yj))
    
    beta1 = yj - yk
    beta2 = yk - yi
    beta3 = yi - yj

    gamma1 = xk - xj 
    gamma2 = xi - xk
    gamma3 = xj - xi

    B = 1/(2*A)*np.array([[beta1, beta2, beta3, 0, 0, 0], [0, 0, 0, gamma1, gamma2, gamma3], [beta1, beta2, beta3, gamma1, gamma2, gamma3]])
    D = E/(1-v**2)*np.array([1, v, 0], [v, 1, 0], [0, 0, (1 - v)/2])

    Ke = t*A*multi_dot([np.transpose(B), D, B])
    return Ke


def Rectangular_Element_SM(E, v, e, IN, Coordinates):
    pass
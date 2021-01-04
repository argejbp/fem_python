import numpy as np
from numpy.linalg import multi_dot
import math

#Solid Mechanics line elements

def Euler_Bernoulli_Truss(E, I, A, xi, yi, xj, yj):

    L = math.sqrt(((xi - xj)**2 + (yi - yj)**2))
    lx = (xi - xj)/L
    ly = (yi - yj)/L
    EI = E*I
    EA = E*A

    ke = np.array([[EA/L, 0, 0, -EA/L, 0, 0], [0, 12*EI/L**3, 6*EI/L**2, 0, -12*EI/L**3, 6*EI/L**2], [0, 6*EI/L**2, 4*EI/L, 0, -6*EI/L**2, 2*EI/L], [-EA/L, 0, 0, EA/L, 0, 0], [0, -12*EI/L**3, -6*EI/L**2, 0, 12*EI/L**3, -6*EI/L**2], [0, 6*EI/L**2, 2*EI/L, 0, -6*EI/L**2, 4*EI/L]])
    MT = np.array([[lx, ly, 0, 0, 0, 0], [-ly, lx, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, lx, ly, 0], [0, 0, 0, -ly, lx, 0], [0, 0, 0, 0, 0, 1]])

    ker = multi_dot([MT.T, ke, MT])
    return ker

def Bar_2D(E, A, xi, yi, xj, yj):

    L = math.sqrt(((xi - xj)**2 + (yi - yj)**2))
    lx = (xi - xj)/L
    ly = (yi - yj)/L
    EA = E*A

    ke = np.array([[EA/L, -EA/L], [-EA/L, EA/L]])
    MT = np.array([[lx, ly, 0, 0], [0, 0, -ly, lx]])

    ker = multi_dot([MT.T, ke, MT])
    return ker

#Heat transfer line elements 

def fin_element(k, h, A, P, xi, yi, xj, yj):
    L = math.sqrt(((xi - xj)**2 + (yi - yj)**2))
    ke = k*A/L*np.array([[1, -1], [-1, 1]]) + h*P*L/6*np.array([[2, 1], [1, 2]])
    return ke

def fin_element_lv(h, To, P, xi, yi, xj, yj):
    L = math.sqrt(((xi - xj)**2 + (yi - yj)**2))
    fe = h*To*P*L/(2)*np.array([[1], [1]])
    return fe 

def conduction_element(k, A, xi, yi, xj, yj):
    L = math.sqrt(((xi - xj)**2 + (yi - yj)**2))
    ke = k*A/L*np.array([[1, -1], [-1, 1]])
    return ke

def convection_element(h, A):
    ke = h*A*np.array([[1, -1], [-1, 1]])
    return ke


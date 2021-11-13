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

def Euler_Bernoulli_3DFrame(E, Izz, Iyy, A, G, J, p1, p2, p3):
    # Coordinates of node i
    xi = p1[0]
    yi = p1[1]
    zi = p1[2]

    # Coordinates of node j
    xj = p2[0]
    yj = p2[1]
    zj = p2[2]

    # Coordinates of reference node
    xk = p3[0]
    yk = p3[1]
    zk = p3[2]

    # Length of the element
    l12 = math.sqrt(((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2))
    # Distance between node i and reference node
    l13 = math.sqrt(((xi - xk)**2 + (yi - yk)**2 + (zi - zk)**2))

    # Main terms of the stiffness matrix
    EA = E*A
    EIzz = E*Izz
    EIyy = E*Iyy
    GJ = G*J

    # Direction cosines in X direction
    l1 = (xj - xi)/l12
    m1 = (yj - yi)/l12
    n1 = (zj - zi)/l12
    Vx = np.array([l1, m1, n1])

    # Direction cosines in Z direction
    V13 = np.array([((xk - xi)/l13), ((yk - yi)/l13), (zk - zi)/l13])
    cross_product = np.cross(Vx, V13)
    Vz = cross_product/np.linalg.norm(cross_product)
    l3 = Vz[0]
    m3 = Vz[1]
    n3 = Vz[2]

    # Direction cosines in Y direction
    Vy = np.cross(Vz, Vx)
    l2 = Vy[0]
    m2 = Vy[1]
    n2 = Vy[2]

    # Transformation matrix
    mt = np.array([[l1,m1,n1],[l2,m2,n2],[l3,m3,n3]])
    MT = np.zeros((12,12))
    MT[0:3,0:3] = MT[0:3,0:3] + mt
    MT[3:6,3:6] = MT[3:6,3:6] + mt
    MT[6:9,6:9] = MT[6:9,6:9] + mt
    MT[9:,9:] = MT[9:,9:] + mt

    # Locations of every individual matrix
    ax = np.array([0,6])
    fzz = np.array([1,5,7,11])
    fyy = np.array([2,4,8,10])
    fxx = np.array([3,9])

    # Adding the terms to the local stiffness matrix
    ke = np.zeros((12,12))
    ke[ax[:,None], ax[:,None]] = ke[ax[:,None], ax[:,None]] + (EA/l12)*np.array([[1, -1], [-1, 1]])
    ke[fzz[:,None], fzz[:,None]] = ke[fzz[:,None], fzz[:,None]] + (EIzz/(l12**3))*np.array([[12,6*l12,-12,6*l12],[6*l12,4*l12**2,-6*l12,2*l12**2],[-12,-6*l12,12,-6*l12],[6*l12,2*l12**2,-6*l12,4*l12**2]])
    ke[fyy[:,None], fyy[:,None]] = ke[fyy[:,None], fyy[:,None]] + (EIyy/(l12**3))*np.array([[12,6*l12,-12,6*l12],[6*l12,4*l12**2,-6*l12,2*l12**2],[-12,-6*l12,12,-6*l12],[6*l12,2*l12**2,-6*l12,4*l12**2]])
    ke[fxx[:,None], fxx[:,None]] = ke[fxx[:,None], fxx[:,None]] + (GJ/l12)*np.array([[1, -1], [-1, 1]])

    # Transforming the local stiffness to global stifness matrix
    ker = multi_dot([np.transpose(MT), ke, MT])
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


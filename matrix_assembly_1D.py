import numpy as np
import stiffness_matrices_1D_Elements as me
def stiffness_matrix_assembly(Problem_Type, coordinates, Geometric_properties, Material_properties, Elements, Ne, Nn):
    #Thermal Problems 
    if Problem_Type == 1:

        KG = np.zeros((Nn, Nn))                                 #Initialiazing Stiffness Matrix
        FG = np.zeros((Nn, 1))                                  #Initialiazing Load Vector    

        for i in range(Ne):

            pos_i = int(Elements[i, 1]) - 1
            pos_j = int(Elements[i, 2]) - 1

            element_dof = np.array([pos_i, pos_j])
            xi = coordinates[pos_i, 0]
            yi = coordinates[pos_i, 1]
            xj = coordinates[pos_j, 0]
            yj = coordinates[pos_j, 1]
            A = Geometric_properties[i, 0]
            
            if Elements[i, 0] == 1:     #Conductive Element
                k = Material_properties[i, 0]
                Ke = me.conduction_element(k, A, xi, yi, xj, yj)    #Calling function to calculate Conductive Element Stiffness Matrix

            elif Elements[i, 0] == 2:   #Convective Element
                h = Material_properties[i, 0]
                Ke = me.convection_element(h, A)                    #Calling function to calculate Convective Element Stiffness Matrix

            elif Elements[i, 0] == 3:   #Fin Element
                p = Geometric_properties[i, 1]
                k = Material_properties[i, 0]
                h = Material_properties[i, 1]
                To = Material_properties[i, 2]
                Ke = me.fin_element(k, h, A, p, xi, yi, xj, yj)         #Calling function to calculate Fin Element Stiffness Matrix
                Fe = me.fin_element_lv(h, To, p, xi, yi, xj, yj)     #Calling function to calculate Fin Element Load Vector
                FG[element_dof, :] = FG[element_dof, :] + Fe          #Adding the element load vector to the global load vector
            else:
                print('Wrong Entry, that element does not exists')
                
            
            KG[element_dof[:, None], element_dof] = KG[element_dof[:, None], element_dof] + Ke    #Adding element stiffness matrix to the global stiffness matrix
    
    #Solid Mechanics Problems
    elif Problem_Type == 2:
        #Coming Soon
        pass

    else:
        print('Wrong Entry, that type of problem does not exists')
import numpy as np

def thermal_boundary_conditions(bc, Elements, KG, FG):
    
    for i in range(len(bc[::,0])):

        tb = bc[i,0]
        value1 = bc[i,2]
        value2 = bc[i,3]

        if tb == 1:     #Fixed temperature (value1 = Temp)
            node = int(bc[i, 1]) - 1
            KG[node, node] = KG[node, node]*1E+10
            FG[node, 0] = KG[node, node]*value1

        elif tb == 2:   #Heat Flow (value1 = Heat Flux; value2 = Area)
            node = int(bc[i, 1]) - 1
            FG[node, 0] = FG[node, 0] + value1*value2

        elif tb == 3:   #Heat convection (value1 = Conv. Coeff; value2 = Temp, value3 = Area)
            node = int(bc[i, 1]) - 1
            value3 = bc[i, 4]

            KG[node, node] = KG[node, node] + value1*value3
            FG[node, 0] = FG[node, 0] + value1*value2*value3

        elif tb == 4:   #Heat generation (value1 = Heat generation; value2 = Lenght of element; value3 = Area)
            element = int(bc[i, 1]) - 1
            nodei = int(Elements[element, 1]) - 1
            nodej = int(Elements[element, 2]) - 1
            value3 = bc[i, 4]

            FG[nodei, 0] = FG[nodei, 0] + value1*value2/2*value3
            FG[nodej, 0] = FG[nodej, 0] + value1*value2/2*value3

        else:
            print("Wrong entry, that boundary condition does not exists")
    
    return KG, FG


def solid_mechanics_boundary_conditions(bc, Elements, KG, FG):

    # for i in range(len(bc[::,0])):
    #     #coming soon
    #     pass

    return KG, FG

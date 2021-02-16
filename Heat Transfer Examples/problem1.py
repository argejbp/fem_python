import matrix_assembly as m_a
import fem_solution as f_s
import boundary_conditions_1d_elements as bc 
import numpy as np

problem_type = 1
Number_of_elements = 2
Number_of_nodes = 3
coordinates = np.array([[0, 0], [0.15, 0], [0.30, 0]])
geometric_properties = np.array([[15], [15]])
material_properties = np.array([[0.9], [0.9]])

bound_cond = np.array([[1, 1, 16, 0, 0], [1, 3, 2, 0, 0]])

elements = np.array([[1, 1, 2], [1, 2, 3]])

KG, FG = m_a.stiffness_matrix_assembly(problem_type, coordinates, geometric_properties, material_properties, elements, Number_of_elements, Number_of_nodes)
KG, FG = bc.thermal_boundary_conditions(bound_cond, elements, KG, FG)
Temps = f_s.steady_state_solution(KG, FG)
print(Temps)

import matrix_assembly as m_a
import fem_solution as f_s
import boundary_conditions_1d_elements as bc 
import numpy as np

problem_type = 1
Number_of_elements = 3
Number_of_nodes = 4
coordinates = np.array([[0, 0], [0.30, 0], [0.45, 0], [0.60, 0]])
geometric_properties = np.array([[1], [1], [1]])
material_properties = np.array([[20], [30], [50]])

bound_cond = np.array([[3, 1, 25, 800, 1], [1, 4, 20, 0, 0]])

elements = np.array([[1, 1, 2], [1, 2, 3], [1, 3, 4]])

KG, FG = m_a.stiffness_matrix_assembly(problem_type, coordinates, geometric_properties, material_properties, elements, Number_of_elements, Number_of_nodes)
KG, FG = bc.thermal_boundary_conditions(bound_cond, elements, KG, FG)
Temps = f_s.steady_state_solution(KG, FG)
print(Temps)

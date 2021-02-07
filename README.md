# fem_python
Finite Element functions in Python for 1D elements (In progress)

These are Python files for solving Finite Element Method problems with 1D Elements (or line elements)
The elements available are: (10/12/2020)
- Heat Transfer:
  - Conduction Element
  - Convection Element
  - Fin Element
  
- Solid Mechanics:
  - Bar element (2d formulation)
  - Euler-Bernoulli Beam (In-plane formulation)
  
It's planned Steady State solutions for Heat Transfer and Solid Mechanics, and Modal Analysis for Solid Mechanics
Changes will be added soon, and it will be updated below this line with date 

02/02/2021 Update

Uploaded files:
  - matrix_assembly_1D.py: It builds up the Stiffness matrix and Load vector (At this time, heat Transfer only)
  - boundary_conditions_1D_elements.py: Modifies the stiffness matrix and load vector to apply boundary conditions (At this time, heat Transfer only)

07/02/2021 Update

Uploaded files:
  - fem_python1D.pdf: It's a file describing the stiffness matrix of every element in matrix_assembly_1D.py
  
A few errors were captured in matrix_assembly_1D.py during tests of Heat Transfer modules, these were modified in the uploaded file


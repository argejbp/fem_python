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

14/02/2021 Update

Uploaded files: 
  - "Heat Transfer example exercises.pdf": It's a file describing two example problemd of 1D Heat Transfer solved by the Finite Element Method, using the files uploaded   here
  - problem1.py: First example explained in "Heat Transfer example exercises.pdf"
  - problem2.py: Second example explained in "Heat Transfer example exercises.pdf"
  
03/07/2021 Update
It's been a while. Here's a new file:
  - "stiffness_matrices_2D_Elements.py": It's a file with functions to calculate stiffness matrix of: 2D Triangle element (for Solid Mechanics problems and Heat Transfer problems) and 3D Axisymmetric Triangle Element (for Heat Transfer)
  
  
  

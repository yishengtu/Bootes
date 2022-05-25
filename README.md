# Bootes
3D Hydrodynamic meshgrid code

Get started
Step 1: write problem generater in src/setup/<problem>.cpp
Step 2: in main.cpp, change #include "setup/<problem>.cpp" from whatever it was to <problem>.cpp
Step 3: in boundary/apply_bc.hpp, include desired boundary conditions
Step 4: "make" under root directory
Step 5: "cd bin"
Step 6: change input.txt to desired inputs
Step 7: "./bootes"


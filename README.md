# Bootes
3D Hydrodynamic meshgrid code

Get started <br>
Step 1: write problem generater in src/setup/\<problem\>.cpp <br>
Step 2: in src/main.cpp, change #include "setup/\<problem\>.cpp" from whatever it was to <problem>.cpp <br>
Step 3: in src/defs.hpp, include only flags needed <br>
Step 4: in src/algorithm/boundary_condition/apply_bc.cpp, include desired boundary conditions <br>
&emsp; Step 4.1: if dust is enabled, also change boundary conditions for dust in boundary_condition/dust/apply_bc_dust.cpp <br> 
Step 5: under root directory, if a directory named "obj" doesn't exist, then create one using "mkdir obj" <br>
Step 6: "make" under root directory <br>
Step 7: "cd bin" <br>
Step 8: change input.txt to desired inputs or use a restart file <br>
Step 9: run the code with <br>
&emsp; Step 9.1: "./bootes.out -i \<input file name\>" if using input file <br>
&emsp; Step 9.2: "./bootes.out -r \<restart file name\>" if using a restart file <br>

Example 1: Kelvin-Homoltz in Cartesian coordinate <br>
Step 1: in src/main.cpp, change #include setup/\<problem\>.cpp" to #include setup/KH.dust.cpp" <br>
Step 2: in src/defs.hpp, enable only Cartesian Coodinate, dust fluid (exclude coagulation) and Density Protection <br>
Step 3: in boundary/apply_bc.cpp, use standard boundary for x1 and x2. x3 can be set to anything (recommand periodic). This step sets boundary condition for hydro <br>
Step 4: in boundary/dust/apply_bd_dust.cpp, use the same boundary conditions as in hydro <br>
Step 5: under root directory, make a directory "obj" <br>
Step 6: "make clean; make" under root directory <br>
Step 7: "cd bin" <br>
Step 8: "./bootes.out -i input.txt.KH" to run the code. <br>

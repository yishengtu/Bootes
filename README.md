# Bootes
3D Hydrodynamic meshgrid code

Get started <br>
Step 1: write problem generater in src/setup/\<problem\>.cpp <br>
Step 2: in main.cpp, change #include "setup/\<problem\>.cpp" from whatever it was to <problem>.cpp <br>
Step 3: in boundary/apply_bc.hpp, include desired boundary conditions <br>
Step 4: under root directory, if a directory named "obj" doesn't exist, then create one using "mkdir obj" <br>
Step 5: "make" under root directory <br>
Step 6: "cd bin" <br>
Step 7: change input.txt to desired inputs or use a restart file <br>
Step 8: run the code with <br>
&emsp Step 8.1: "./bootes.out -i \<input file name\>" if using input file <br>
&emsp Step 8.2: "./bootes.out -r \<restart file name\>" if using a restart file <br>


# AssemblyPlanFromObservation
Implementation of Core of Assembly Plan from Observation (APO)

[1] J. TaKamatsu, K. Ogawara, H. Kimura, and K. Ikeuchi, "Recognizing Assembly Tasks Through Human Demonstration", IJRR, vol. 26, no. 7, pp.641-659, 2007.

Currently, I implemented to 
- formulate the feasible infinitesimul displacement from the contact information (type, contact position, outer normal, etc.)
- calculate the index of infinitesimul displacement (maintaining, detaching, constraining, restricted DOF in tanslation/rotation)

## Requirement
- CMake
- Polyhedral Convex Cones (https://github.com/j-taka/PolyhedralConvexCones) and its dependencies

## Compile
1. Download all requirements
2. Compile Polyhedral Convex Cones
3. Change PCC_PATH in CMakeList.txt correctly
4. Use cmake to compile the program

## Test
just run APO_Recognition_sample.exe

## Author
Jun Takamatsu j-taka@is.naist.jp

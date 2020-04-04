# AssemblyPlanFromObservation
Implementation of Core of Assembly Plan from Observation (APO)

[1] J. TaKamatsu, K. Ogawara, H. Kimura, and K. Ikeuchi, "Recognizing Assembly Tasks Through Human Demonstration", IJRR, vol. 26, no. 7, pp.641-659, 2007.

Currently, I implemented to 
- formulate the feasible infinitesimul displacement from the contact information (type, contact position, outer normal, etc.)
- calculate the index of infinitesimul displacement (maintaining, detaching, constraining, restricted DOF in tanslation/rotation)
- correct vision errors using roughly estimated contact states
- generate a trajectory from contact state transitions.

## Requirement
- CMake
- [Polyhedral Convex Cones](https://github.com/j-taka/PolyhedralConvexCones) and its dependencies
- [OpenCascade](https://www.opencascade.com/)
- [Qt](https://www.qt.io/)
- Boost

## Compile
1. Download all requirements
2. Compile Polyhedral Convex Cones
3. Change file passes in CMakeList.txt correctly
4. Use cmake to compile the program

## Test
- run APO_Recognition_sample.exe (calculate motion DOFs)
- run VisionErrorCorrection.exe error_correct.pose.list peg_in_hole.pose.list (vision error correction) and open error_correct_pose.list by [APOViewer](https://github.com/j-taka/APOViewer)
  - Necessary files (peg.step, hole.step, peg_in_hole.pose.list) are included in APOViewer
- run TaskAnalyzer.exe state.pose.list error_correct.pose.list. You can see the list of the contact states and motion DOFs on the terminal and see the actual situations of the states using [APOViewer](https://github.com/j-taka/APOViewer)
- run MotionGenerationFromStateTransitions.exe motion.pose.list state.pose.list init.pose.list and open motion.pose.list by APOViewer

## Author
Jun Takamatsu j-taka@is.naist.jp

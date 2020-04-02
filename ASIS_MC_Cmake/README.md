This part is the work on importance sampling + adaptive stratification.
1. Importance Sampling: Gradient Based Optimization for finding the optimal change of measure.
2. Adaptive Stratification:
3. Apply to toy example:Normal(0,1) and Asian Option Pricing.


##### Important: How to compile and run?####
0. Before compiling, you should find 4 folders: "include"---containing the header files; "src"---containing the cpp files; "bin"---empty folder where we shall store the executable after compiling; "build"---empty folder where we shall store the make file after we run "cmake"
1. Enter the empty "./build" directory by "cd ./build"  (if ./build is not empty, delete all files in it)
2. press "cmake .." to generate the make file
3. press "make" to compile codes in "../src"
4. The runnable files are in "../bins"
5. press "../bins/asianOptionNew" to run the asian option pricing experiments and the output of estimating would be stored in "call_res.csv" and "put_res.csv"
6. press "../bins/simpleNormal" to run the toy example described in our report. 
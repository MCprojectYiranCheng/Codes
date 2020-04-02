##This part is the work on importance sampling + adaptive stratification.

1. Importance Sampling: Gradient Based Optimization for finding the optimal change of measure.
2. Adaptive Stratification:
3. Apply to toy example:Normal(0,1) and Asian Option Pricing.


##Important: How to compile and run?
0. Before compiling, you should find 2 folders: "include"---containing the header files; "src"---containing the cpp files; and also 1 file "CMakeLists.txt" helping generating the make file.

1. press "cmake ." to generate the make file

2. press "make" to compile codes in "./src"

3. The runnable files are in "./bin"

4. press "./bin/asianOptionNew" to run the asian option pricing experiments and the output of estimating would be stored in "call_res.csv" and "put_res.csv"

5. press "./bin/simpleNormal" to run the toy example described in our report. 
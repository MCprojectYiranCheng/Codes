This part is the work on randomized QMC. To do our job, we have created mainly three classes including shifted_halton, Gaussian_SH and RQMC_Gaussian. The function for these three classes are as follows:
1. shifted_halton : to generate the shfited halton sequence by providing the dimension and shifted value.
2. Gaussian_SH    : to generate a sequence to simulate a sequence of i.i.d Gaussian variable by using shited halton sequence. We have used the Box-Muller method to generate guassian.
3. RQMC_Gaussian  : to do Randomized QMC for calculating $\mathbb{E[f(X]]$ where X follows a standard d-dimension gaussian distribution. 



 We have two tasks for this part.
1. Write a class for generate randomized sequence for RQMC
2. Do the simulation for N(0,1) and discretized option pricing


# 1. Step by step guide to install our tool for doing Randomized QMC

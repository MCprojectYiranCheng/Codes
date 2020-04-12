## Bref introuduction

This part is the work on randomized QMC. To do our job, we have created mainly three classes including shifted_halton, Gaussian_SH and RQMC_Gaussian. The function for these three classes are as follows:
1. shifted_halton : to generate the shfited halton sequence by providing the dimension and shifted value.
2. Gaussian_SH    : to generate a sequence to simulate a sequence of i.i.d Gaussian variable by using shited halton sequence. We have used the Box-Muller method to generate guassian.
3. RQMC_Gaussian  : to do Randomized QMC for calculating E[f(X)] where X follows a standard d-dimension gaussian distribution. Return the estimated result and the half length of the 95% confidence interval.

In the main.cpp, we have wrote the code to use our class to generate the data for doing a randomized QMC for estimating a standard gaussian varible's mean and for pricing the arithmemaic Asian call and put option.

## A. Step by step guide
Run the following bash code in the terminal under the directory of `Codes/RQMC/` in order to run our code for doing randomized QMC job.
```Bash
mdkir build
cd build
cmake ..
make
../bin/RQMC
```

## B. How to use our class RQMC_Gaussian?
RQMC_Gaussian is used to do Randomized QMC for solving the problem in the for $ \sum_{\forall i}{x_i^{2}} $

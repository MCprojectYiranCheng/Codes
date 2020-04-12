## Bref introuduction

This part is the work on randomized QMC. To do our job, we have created mainly three classes including shifted_halton, Gaussian_SH and RQMC_Gaussian. The function for these three classes are as follows:
1. shifted_halton : to generate the shfited halton sequence by providing the dimension and shifted value.
2. Gaussian_SH    : to generate a sequence to simulate a sequence of i.i.d Gaussian variable by using shited halton sequence. We have used the Box-Muller method to generate guassian.
3. RQMC_Gaussian  : to do Randomized QMC for calculating E[f(X)] where X follows a standard d-dimension gaussian distribution. Return the estimated result and the half length of the 95% confidence interval.

In the main.cpp, we have wrote the code to use our class to generate the data for doing a randomized QMC for estimating a standard gaussian varible's mean and for pricing the arithmemaic Asian call and put option.

## A. Step by step guide
Run the following bash code in the terminal under the directory of `Codes/RQMC/` in order to run our code for doing randomized QMC job.
```Bash
mkdir data
mdkir build
cd build
cmake ..
make
../bin/RQMC
```

## B. How to use our class RQMC_Gaussian?
RQMC_Gaussian is used to do Randomized QMC for solving the problem in the form of E(f(X)) where X is a d-dimension random vertor whose law is N(0,I_d) where I_d is a d-dimension identity matrix. So in order to use our class, and f is a fonction from R^d -> R. So we need to know the defintion of f in our specific problem in order to use RQMC_Gaussian.
1. The type of f function needs to be std::function<double(std::vector<double>)>
2. N,I are the length of low discrepancy sequence and the number of copied of shifted sequence.
3. d is the dimension of problem
4. gen is a random generator
Then we can declare an object of RQMC_Gaussian class by `RQMC_Gaussian object(I,N,d,f,gen)`. Once we execute `object()`, it will do the randomized QMC job. Then the results are saved in the memeber of object, `object.estimator` and `object.half_CI`. We can use the following `to_csv` function to save the results into a csv file.


``` C++
void to_csv(std::string dir, const RQMC_Gaussian rqmc){
  std::ofstream myfile;
  myfile.open (dir);
  myfile << rqmc <<"\n";
  myfile.close();
}
```

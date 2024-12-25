# scFMA

`scFMA` (single-cell flexible methylation analysis) is a statistical method that can be used to analyze DNA methylation patterns in 
single-cell bisulfite sequencing data. In this method, discrete approximation technique and MM algorithm are introduced, which can 
construct wald statistic based on the binomial model with flexible random effects to automatically identify the differential methylation 
regions (DMRs) between cells. Thus, `scFMA` can avoid the complicated procedure of model estimation, and even solve the problem of estimation 
bias caused by misspecifying the distribution of random effect in advance. 

## Install R package
Here is a R package called `scFMA.0.0.0.9000` that you need to download and install the package in R. After librarying the 'scFMA' package, you can call 
"sim_scFMA" function and "real_scFMA" function to perform differential methylation analysis on the two types of samples in simulated and real data, 
respectively. The output of "sim_scFMA" function are p-values, regional level methylation difference (Δ) between and the parameter estimation results 
of the methylation distribution between the two types of samples. And the output of "real_scFMA" function is p-values and regional level methylation 
difference (Δ). By setting different combinations of thresholds, the DMRs that can be finalized based on the outputs of our method should satisfy the 
p-values not exceeding the p-value cutoff and ∆ exceeding the ∆ cutoff. 

## Simulation study
In simulation experiments, our approach considers differential methylation analysis under two random effects distributions to illustrate the 
broad applicability of `scFMA`. 

### Zero-one inflated beta distribution
Run `generate4c8c_beta.r` and `generate8c8c_beta.r` to simulate 0-1 inflated beta-binomial model on human-generated data and estimate the model. A list 
called `testRegion` consisting of methylation reads and total reads for the two types of samples, the parameter estimates `res.8c` and `res.4c` for the 
two types of samples, and the random numbers `pio_list` generated by the 0-1 inflated beta distribution are then available. 
```
### difference experiment
run: generate4c8c_beta.r
### indifference experiment
run: generate8c8c_beta.r
```

### Zero-one inflated simplex distribution
Run `generate4c8c_simplex.r` and `generate8c8c_simplex.r` to simulate 0-1 inflated simplex-binomial model on human-generated data and estimate the model. 
```
### difference experiment
run: generate4c8c_simplex.r
### indifference experiment
run: generate4c8c_simplex.r
```
### Model estimation
Subsequently, based on the simulation data generated from different cases, runing the `sim_scFMA` function allows for the determination of DMRs under 
different threshold scenarios, as well as the visualisation of the estimates of the distribution of random effects.
```
load("difference.beta_est.Rdata") #### nodifference.beta_est.Rdata  #difference.simplex_est.Rdata  #nodifference.simplex.Rdata
testRegion <-res$testRegion
res.8c <- res$res.8c
res.4c <- res$res.4c
pio_list <- res$pio
sim_result <- sim_scFMA(testRegion,res.8c, res.4c,pio_list,n1=48,n2=25,start=1,end=10)
pvalue <- sim_result[1,]
delta <- sim_result[2,]
pai8c <- sim_result[4:33,]
pai4c <- sim_result[35:59,]
```
## Data analysis
Run the `real_scFMA` function to identify DMRs in regions filtered out of the real data.
```
load("chr5_4cVS8c.Rdata")
treadn <- chr5_data$treadn
treadx <- chr5_data$treadx
testRegion <- chr5_data$testRegion
sample8c <- c(2:12,16,17,19,21:25,27,29,33,35,36,38:40,43:45,47:49,52:55,57:60,62:64,66,69,70,73)
sample4c <- c(1,13:15,18,20,26,28,30:32,34,37,41,42,46,50,51,56,61,65,67,68,71,72)
real_result <- real_scFMA(testRegion,treadn, treadx,sample8c,sample4c,n1=48,n2=25)
pvalue <- real_result[1,]
delta <- real_result[2,]
```
We put some of the data and code needed for the simulation experiments are placed in the simultaion folder. And the raw data and filtering procedure 
for the real data are placed in the realdata folder.

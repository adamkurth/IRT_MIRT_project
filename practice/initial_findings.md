## Initial Findings

Assuming \(n = 500\).

### Directions: 

In R's `mirt` package: 
-  **2PL Model Use**: simulate data from a 2-parameter logistic (2PL) model from Item Response Theory (IRT).
-  **True Parameter Specificiations**: true parameters for item characteristic curves (\(a = \text{ discrimination }, b = \text{ item difficulty }\)) are given. For the 2PL model, these parameters directly influence the shape and position of these curves. 
- **Skewed and Normal Distributions**: simulate response data based on the left-skewed, right-skewed, and standard normal distributions.
- **Crossing Parameters for 16 Items**: Creating combinations of `a` and `b` parameters for 16 items, so that each unique pair is treated as a seperate "item" in the analysis. Since there are 4 values for each parameter, this results in 16 unique combinations.
- **Reporting Requirements**: 
  - **RMSE**: Root Mean Squared Error (RMSE) for each item, averaged over 100 replications.
  - **Bias**: Bias for each item, averaged over 100 replications. 
  - **Coverage Time**: Time to converge for each item.

### Step 1:

1. *Generate Skewed and Normal Distributions*: For each distribution type (left, right and standard normal), simulate response data based on the 2PL model with given `a` and `b` parameters.
2. *Crossing `a`, `b` Parameters*: Create a dataset for each of the 16 items by crossing `a` and `b` parameters. Each will have a unique combination of `a` and `b` parameters and then simulate responses for these items under each distribution type.
3. *Replications:* For accurate estimation, conduct 100 replications of data simulation and parameter estimation for each item under each distribution type. This involves simulating new datasets and estimating parameters 100 times to compute the average RMSE, bias, and coverage time for each item.
4. *Parameter Estimation*: Utilize the `mirt` package to estimate the `a` and `b` parameters from the simulated response data. 
5. *Compute RMSE and Bias*: After parameter estimation, for each replicaiton, compute the RMSE and bias for each item. Average these values over the 100 replications. 
6. *Computation Time Tracking*: Record the computation time for each replication to report the average time taken for parameter estimation across replications. 

### Step 2:

- Aggregate the results across all items, distribution types, and replications. Compute the average RMSE, bias, computation time into a structured format. 
- Use multidimensional arrays or lists to organize results before summarizing them into readable format.
- Present the results in a table or report that clearly distinguishes the different estimation methods and density types, and distribution types used.

#### Psuedo Code:

``` r
# simulate_data is function that outputs response_data
# estimate_parameters is function that outputs estimatedd parameters and computation time.
results <- list()
for(distribution in c("left.skew", "right.skew", "stnd.norm")){
  for(a_val in c(0.5, 1, 1.5, 2)){
    for(b_val in c(0, -1, 1, 0)){
      simulate_data <- simulate_data(distribution, a_val, b_val)
      est.param, comp.time <- estimate_parameters(simulate_data)
      # Compute RMSE and Bias and store (with comp.time) in results
      results[[paste(distribution, a_val, b_val)]] <- list(est.param, comp.time)
    }
  }
}

```


### 1. Examine the following distributions:
  - Skew Right
  - Skew Left
  - Standard Normal

### 2. Examine Estimation Methods:
 - **Assuming all other arguments are default.**
 
 - *Block & Lieberman approach (BL)*:
    - From R documentation: "Block and Lieberman (1974) approach estimates the latent distribution using a nonparametric kernel density estimator. This method is only applicable for unidimensional models estimated with the EM algorithm. By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm"
    - Call:
       ``` r
       mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='BL',dentype='Gaussian')
       # Converged within 1e-08 tolerance after 99 BL iterations.
       ```

 - *Expectation- Maximization (EM) algorithm: (default)*
   - From R Documentation:"The 'EM' method estimates the latent distribution using the Expectation-Maximization (EM) algorithm described by Bock and Aitkin (1981). This method is only applicable for unidimensional models estimated with the EM algorithm. By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm"
   - Call:
       ``` r
       mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='EM', dentype='Gaussian')
       #FAILED TO CONVERGE within 1e-04 tolerance after 500 EM iterations.
       ```
 - *Metropolis-Hastings Robbins-Monro (MHRM):*
   - From R documentation: "The 'MHRM' method estimates the latent distribution using a Metropolis-Hastings Robbins-Monro (MHRM) algorithm described by Cai and Hansen (2013). This method is only applicable for unidimensional models estimated with the EM algorithm. By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm"
   - Call:
     ``` r
     mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='MHRM',dentype='Gaussian')
     # Converged within 0.001 tolerance after 75 MHRM iterations.
     ```


- *Monte Carlo EM (MCEM):* 
  - From R Documentation: "The 'MCEM' method estimates the latent distribution using a Monte Carlo EM (MCEM) algorithm described by Cai and Hansen (2013). This method is only applicable for unidimensional models estimated with the EM algorithm. By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm"
  - Call:
    ``` r
    mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='MCEM',dentype='Gaussian')
    # FAILED TO CONVERGE within 1e-04 tolerance after 500 MCEM iterations.
    ```

- *Stochastic EM (SEM):*
  - From R Documentation: "The 'Stochastic EM' method estimates the latent distribution using a stochastic EM algorithm described by Cai and Hansen (2013). This method is only applicable for unidimensional models estimated with the EM algorithm. By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm"
  - Call:
    ``` r
    mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='SEM',dentype='Gaussian')
    # Converged within NA tolerance after 100 SEM iterations.
    ```

- *Quasi-Monte Carlo EM (QMCEM):*
  - From R documentaiton: "The 'QMC-EM' method estimates the latent distribution using a quasi-Monte Carlo EM algorithm described by Cai and Hansen (2013). This method is only applicable for unidimensional models estimated with the EM algorithm. By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm"
  - Call:
    ``` r
    mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='QMCEM',dentype='Gaussian')    
    # FAILED TO CONVERGE within 1e-04 tolerance after 500 QMCEM iterations.
    ```


### 3. Examine Density Type:
 - *Gaussian* (default):
   -  From R documentation: "'Gaussian' (default) assumes a multivariate Gaussian distribution with an associated mean vector and variance-covariance matrix"
   - Call:
        ``` r
        mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='EM',dentype='Gaussian')
        # FAILED TO CONVERGE within 1e-04 tolerance after 500 EM iterations.
        ```

 - *EH*:
   - From R documentation: "'empiricalhist' or 'EH' estimates latent distribution using an empirical histogram described by Bock and Aitkin (1981). Only applicable for unidimensional models estimated with the EM algorithm. For this option, the number of cycles, TOL, and quadpts are adjusted accommodate for less precision during estimation (namely: TOL = 3e-5, NCYCLES = 2000, quadpts = 121)"
   - *"Error: Too few degrees of freedom. There are only 15 degrees of freedom but 128 parameters were freely estimated."*
   - Call:
        ```r
        mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='EM',dentype='EH')
        # Error: Too few degrees of freedom. There are only 15 degrees of freedom but 128 parameters were freely estimated.
        ```

 - *EHW*: 
   - From R documentation: "empiricalhist_Woods' or 'EHW' estimates latent distribution using an empirical histogram described by Bock and Aitkin (1981), with the same specifications as in dentype = 'empiricalhist', but with the extrapolation-interpolation method described by Woods (2007). NOTE: to improve stability in the presence of extreme response styles (i.e., all highest or lowest in each item) the technical option zeroExtreme = TRUE may be required to down-weight the contribution of these problematic patterns"
   - *"Error: Too few degrees of freedom. There are only 15 degrees of freedom but 126 parameters were freely estimated."*
   - Call:
        ```r
        mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='EM',dentype='EHW')
        # Error: Too few degrees of freedom. There are only 15 degrees of freedom but 126 parameters were freely estimated.
        ```

  - *Davidian-#*
    - From R documentation: "'Davidian-#' estimates semi-parametric Davidian curves described by Woods and Lin (2009), where the # placeholder represents the number of Davidian parameters to estimate (e.g., 'Davidian-6' will estimate 6 smoothing parameters). By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm
    - Call: "Davidian-2", "Davidian-4", "Davidian-6".
        ```r
        mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='EM',dentype='Davidian-2')
        # FAILED TO CONVERGE within 1e-04 tolerance after 500 EM iterations.

        mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='EM',dentype='Davidian-4')
        # FAILED TO CONVERGE within 1e-04 tolerance after 500 EM iterations.

        mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='EM',dentype='Davidian-6')
        # FAILED TO CONVERGE within 1e-04 tolerance after 500 EM iterations.
        ```

4. Outcome Metrics: 
   - RMSE (for each item, averaged over 100 replications)
   - Bias (for each item, averaged over 100 replications)
   - Coverage time


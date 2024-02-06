## Initial Findings

Assuming \(n = 500\).

### Directions: 

1. Examine the following distributions:
  - Skew Right
  - Skew Left
  - Standard Normal

2. Examine Estimation Methods:
 - **Assuming all other arguments are default.**
 
 - *Block & Lieberman approach (BL)*:
    - From R documentation: "Block and Lieberman (1974) approach estimates the latent distribution using a nonparametric kernel density estimator. This method is only applicable for unidimensional models estimated with the EM algorithm. By default, the number of quadpts is increased to 121, and this method is only applicable for unidimensional models estimated with the EM algorithm"
  - Call:
     ``` r
     mirt(data=response.dataframes$left.skew, model=1, itemtype='2PL', method='BL',dentype='Gaussian')
     # Converged within 1e-08 tolerance after 99 BL iterations.
     ```

 - *Expectation- Maximization (EM) algorithm: (default)*
   - From R Documentation: "The default is 'EM', for the standard EM algorithm"
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


3. Examine Density Type:
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


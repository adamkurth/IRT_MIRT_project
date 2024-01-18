# What Am I Doing?! 

This is just a quick summarization to understand what the MIRT package in R actually does, and some things that I should focus on to get started. 

# First, 
## Overarching Understanding

- **Objective**: The goal is to estimate the parameters (a, b, c) of the test items based on the responses of the individuals. These parameters are then used to understand how different items function in measuring latent traits (ability/aptitude).

- **Data Structure**: Your dataset where each row is an individual and each column is an item (question on test).

- **Model**: Model fitting involves finding the best set of item parameters that explain the observed responses, this is where we optimize. 

# Second,
## Optimization in MIRT

- Optimizing is at the heart of MIRT model, the idea is to find the parameters that maximize the likelihood of observing the data.

- **Log-Likelihood**: Statistical measure used to assess how well the model (with its given parameters) explains the data. The aim is to maximize the log-likelihood.

- **Algorithm**: Algorithms like the Newton-Raphson algorithm, the EM algorithm is used for optimization. These algorithms iteratively adjust parameters to find max log-likelihood. Each iteration involves calculating the log-likelihood and adjusting the parameters (to maximize) based on the algorithm's rules. 

# Third,
## Key Arguments in `mirt` 

- **method**: This argument type specifies the method to be used. Common methods include 'MH-RM' (Metropolis-Hastings Robbins-Monro), 'EM' (Expectation-Maximization), etc.

- **optimizer**: This argument type specifies the numerical iterative optimizer to be used, such as Newton-Raphson, or 'BDGS', etc. These algorithms have different ways of navigating the parameter space to find the max log-likelihood.

- **dentype**: This argument refers to the density type used in the estimation, especially relevant for certain estimation methods, such as 'Gaussian' (multivariate Gaussian dist), 'empiricalhist' (empirical histogram), 'Davidian-#' (semi-parametric Davidian curves), etc. It affects how the log-likelihood is calculated, and can influence the optimization process.

# Summary
When fitting a MIRT model, you're engaging in an optimization process where you're looking to find the best set of item parameters to explain the data. This involves choosing appropriate estimation methods, an optimizer for finding the max log-likelihood, and a density type for calculating likelihoods. Each choice affects the way the model navigates the complex parameter space to arrive at the most accurate, and reliable estimates of item parameters.

# MIRT Iterations

This repository contains scripts for running simulations using Multidimensional Item Response Theory (MIRT). The scripts are written in R and are used to generate and analyze response data.

## Goals

The main goals of this repository are:

1. Generate response data based on different distributions (standard normal, left skew, right skew). 
2. Simulate response data and write it to .txt files.
3. Replicate this analysis for different software packages (flexMIRT, mirt).
4. Read the response data from the .txt files and analyze it using R's `mirt`, or `flexMIRT` packages. 

## File Structure
.
|__ .gitignore
|__ flex/
|   |__ examples/
|   |__ practice_flex/
|   |__ scripts/
|__ mirt/
|   |__ demo/
|   |__ findings/
|   |__ iterations/
|   |__ markdown/
|__ parallel_fit_mirt.log
|__ practice/
|   |__ 2PL_demo/
|   |__ 3PL_demo/
|   |__ alter_theta_dist.r
|   |__ findings/
|   |__ initial_findings.md
|   |__ iterations/
|   |__ markdown/
|   |__ other/
|   |__ revised_sim_2.r
|   |__ test.r
|   |__ theta_skew.r
|__ README.md
|__ response_data/
    |__ response_left_skew_rep_01.txt
    |__ response_left_skew_rep_02.txt
    |__ response_left_skew_rep_03.txt
    ...
    |__ response_right_skew_rep_01.txt
    |__ response_right_skew_rep_02.txt
    |__ response_right_skew_rep_03.txt
    ...
    |__ response_standard_normal_rep_01.txt
    |__ response_standard_normal_rep_02.txt
    |__ response_standard_normal_rep_03.txt
    ...

#### **flex/**
The `flex/` directory contains several subdirectories, each with a specific purpose:

- `flex/scripts/`: Contains scripts for running MIRT analysis using the `flexMIRT` package. If you're looking to run MIRT analysis using `flexMIRT`, this is where you'd start. 

- `flex/scripts/run.r`: This script runs the MIRT analysis using the `flexMIRT` package. `flexMIRT` must be installed and proper liscence, `run.r` only accesses the script given and runs the analysis.

- `flex/examples/`: This directory contains examples of using the flexMIRT package for MIRT analysis. It's a good place to start if you're new to flexMIRT. This contains all downloaded documents available from flexMIRT website.

- `flex/practice_flex/`: Contains practice scripts for using the `flexMIRT` package. If you're new to `flexMIRT`, this is a good place to understand the new ideas if you're new to flexMIRT.

#### **mirt/**
The mirt/ directory contains several subdirectories, each with a specific purpose:

- `mirt/findings/`: This directory contains the findings and results from the MIRT simulations. It's where you'll find the output data from your simulations.

- `mirt/iterations/`: This directory contains scripts for running MIRT simulations. If you're looking to modify the simulation parameters or the simulation itself, this is where you'd do it.

- `mirt/markdown/`: This directory contains markdown files with notes and summaries of MIRT concepts. It's a good place to look if you're trying to understand the underlying theory behind the simulations.

- `mirt/parallel_fit_mirt.log`: This is a log file for parallel fitting of MIRT models. It's useful for debugging and understanding the performance of the parallel fitting process.

- `mirt/response_data/`: This directory contains response data generated from the MIRT simulations. It's where the raw output of each simulation is stored.

- `mirt/demo/`: This directory contains demo scripts for running MIRT analysis using the mirt package. If you're new to MIRT or this codebase, running these demos is a good first step.


#### **mirt/iterations/**
- `mirt/iterations/iter_4_sim.r`: Contains the main simulation script. It includes functions for generating response data, simulating responses, and exporting the data to .txt files.

- `mirt/iterations/iter_5_sim.r`: An updated version of the simulation script with additional functionality.

- `mirt/iterations/iter_6_sim.r`: The latest version of the simulation script with the most recent updates and improvements.

- `mirt/iterations/tempCodeRunnerFile.r`: A temporary file used for running code snippets (in Visual Studio Code).

- `practice/iterations/iter_4_sim.r`: A practice version of the simulation script used for testing and development.

#### **mirt/demo/**
The `mirt/demo` directory contains R scripts demonstrating the usage of the Multidimensional Item Response Theory (MIRT) package. Here's a brief description of the files:

- `explain.r`: This script explains the parameters of the 2PL (2 Parameter Logistic) model, which includes discrimination and difficulty parameters. It uses the `mirt::coef()` function to extract the estimated parameters from the MIRT analysis. The output includes a dataframe with estimated item parameters and vectors that contain the estimated mean and covariance matrix of the latent trait.

- `2PL_demo.r`: This script records the time to convergence and estimation accuracy of the 2PL model. It uses the `mirt::mirt()` function to fit a maximum likelihood factor analysis model to the data and the `mirt::coef()` function to extract the estimated parameters from the MIRT analysis.

## How to Run

To run the simulations, you can use the `source` function in R to load the script, and then call the `export.data` function with the desired parameters. For example:

```r
# latest version of the simulation script is 'iter_6_sim.r'
source('mirt/iterations/iter_6_sim.r')
export.data(cross.param, output.dir='response_data', n=300, replications=10, seed=123)
```

This will generate response data, simulate responses, and write the data to .txt files in the `response_data` directory.
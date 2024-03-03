# IRT Simulations Project

This project contains a collection of R scripts and data files for simulating and analyzing response data based on different distribution types.

## R's Multiple Item Response (MIRT) Model:
### Directory Structure

- `practice/iterations/`: Contains R scripts for simulating response data and performing analysis.
- `flex/examples/`: Contains example data files.

### Key Files

- `iter_4_sim.r`: This script generates response data based on different distribution types, writes the data to .txt files, and performs analysis using the `mirt` package.
- `iter_5_sim.r`: This script extends the functionality of `iter_4_sim.r` by adding more distribution types and refining the data simulation process.
- `fitRHO-DINA_Interaction-irt.txt`: An example data file.
- `fitRHO-DINA_Interaction-ssc.txt`: Another example data file.

### Key Functions

- `quick.gen.dist`: Generates theta values based on a specified distribution type.
- `quick.sim.response`: Simulates response data based on theta values and cross parameters.
- `read.data`: Reads a specific response data file.
- `export.data`: Simulates response data for different distribution types and writes the data to .txt files.
- `fit.mirt`: Performs analysis on the response data using the `mirt` package.

### How to Run

To run the scripts, open them in your R environment and execute them. Make sure to set your working directory to the project root.

### Dependencies

This project requires the `mirt` R package. Install it with `install.packages("mirt")`.

 ---

## FlexMIRT Model:

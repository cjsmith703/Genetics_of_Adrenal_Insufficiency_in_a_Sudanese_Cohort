# Genomics Interpolator

## Overview

**Genomics Interpolator** is a command-line Python tool designed to interpolate missing values in a matrix using the mean of non-diagonal neighboring values. This tool requires [Conda](https://docs.conda.io/en/latest/) to be pre-installed.

## Download and Installation

After downloading the `genomics_interpolator` directory, follow the steps below to set up and install the tool:

### 1. Navigate to the project directory

```bash
# Move to the genomics_interpolator directory
cd genomics_interpolator
```
### 2. Set up the Conda environment 

```bash
# create a conda environemt with the required dependancies using the environment.yml file
conda env create -f environment.yml
conda activate genomics-env
```

This will install the following dependencies:

python=3.10

pandas

numpy

pyfiglet

### 3. Install Genommics Interpolator 

```bash
# install genomics-interpolator to run on the command line
pip install . 
```

## Usage 
The directory includes a data subdirectory that contains an example CSV file. To run the Genomics Interpolator on the command line, use the following command:

```bash
# run on the command line
genomics-interpolator data\input_test_data.csv
```

The output file will be saved in the same directory as the input file, with _interpolated appended to the filename. For example:
```bash
data/input_test_data.csv → data/input_test_data_interpolated.csv
```
# Fast-GPU-PCC

This repository contains the implementation of Fast-GPU-PCC algorithm.


# Compilation:
Please use the following command for compiling the code:

nvcc -lcublas -O2 -arch=compute_35 -code=sm_35 -std=c++11 CPU_side.cpp GPU_side.cu -o out


# Running

Place the time series data in a text file called data.txt in the current directory. This file will be the input of the Fast-GPU-PCC method.

To run the code use the following command:

./out N M

M: number of voxels

N: length of time series

The resulting correlations will be stored in a binary file called corrs.bin.

# Note

- The data should be stored in row major format (first N elements corresponds to time series of first element, second N elements corresponds to time series of first element and …)

- Example of an input file is having 5 voxels with time series of length 3 is shown in data.txt

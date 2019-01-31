# Fast-GPU-PCC

This repository contains Fast-GPU-PCC algorithm


# Compilation:

nvcc -lcublas -O2 -arch=compute_35 -code=sm_35 -std=c++11 CPU_side.cpp GPU_side.cu -o out


# Running
./out N M 

M: number of elements

N: length of time series

#note:

Make sure that you add the path to your data in the source code

The data should be stored in row major format (first N elements corresponds to time series of first element, second N elements corresponds to time series of first element and …)


Example of an input file is having 5 voxels with time series of length 3 is shown in example_data.txt

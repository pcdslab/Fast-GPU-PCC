# Fast-GPU-PCC

This repository contains the implementation of  **Fast-GPU-PCC** algorithm. This algorithm is used for computing pairwise Pearson’s correlation coefficient. Based on the symmetric property of Pearson’s correlation, this approach returns N(N-1)/2 correlation coefficients located at strictly upper triangle part of the correlation matrix.

Please use the following instructions for compiling and running the code.

## Compilation

Use the following command to compile the code.

```
nvcc -lcublas -O2 -arch=compute_35 -code=sm_35 -std=c++11 CPU_side.cpp GPU_side.cu -o exec
```

## Running

Place the fMRI time series data in a text file called data.txt in the current directory. This file will be the input of the Fast-GPU-PCC method.

To run the code use the following command:

```
./exec N M W
```

|Argument |Description|
|--|--|
|**N** |Number of voxels|
|**M** |Length of time series|
|**W** |Option for writing the results into the binary or text file, possible values: {b: binary, t: text}|


The correlations will be stored in a binary file called corrs.bin or a text file called corrs.txt based on the value of parameter W.

## Note

The data should be stored in the row major format (first N elements corresponds to the time series of first element, second N elements corresponds to the time series of second element and …)


Example of an input file containing 5 voxels with time series of length 3 is shown in the file data.txt.


## Usage example
Use the following commands for computing pariwise corrlations of the sample file (data.txt) and storing the result into a binary file:

```
nvcc -lcublas -O2 -arch=compute_35 -code=sm_35 -std=c++11 CPU_side.cpp GPU_side.cu -o exec

./exec 5 3 b
```

As the result of running this example, 10 pairwise correlation (upper triangle part of the correlation matrix) will be stored in the binary file called corrs.bin.

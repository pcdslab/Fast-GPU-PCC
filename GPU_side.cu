#include "cublas_v2.h"
#include "cuda_runtime.h"
#include "memory.h"
#include <iostream>
#include <ctime>
#include<stdio.h>
#include <string.h>
#include <iomanip>
#include <fstream>
#include <stack>
#include<sstream>
#include<math.h>
using namespace std;

long long remaining_N2(int , int ,long long );
long long remaining_N(int , int ,int );
void preprocessing(float * , int ,int );
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

__global__ void ker(float * cormat, float * upper,int n1,int n)
{
long idx = blockDim.x*blockIdx.x+threadIdx.x;
long i = idx%n1;
long j = idx/n1;
if(i<j && i<n1 && j<n)
{
        long tmp=i;
        tmp*=(i+1);
        tmp/=2;
        long tmp_2=i;
        tmp_2*=n;
        tmp_2=tmp_2-tmp;
       tmp_2+=j;
       tmp_2-=i;


upper[tmp_2-1]=cormat[j*n+i];
}
}


__global__ void ker2(float * cormat, float * upper,int n1,int n,long long upper_size,int N,int i_so_far,long long M1)
{
long long idx = blockDim.x;
idx*=blockIdx.x;
idx+=threadIdx.x;
long i = idx/n;
long j = idx%n;

if(i<j && i<n1 && j<n)// &&i<N &&j<N && idx<(n1*n))
{
        long long tmp=i;
        tmp*=(i+1);
        tmp/=2;
        long long tmp_2=i;
        tmp_2*=n;
        tmp_2=tmp_2-tmp;
        tmp_2+=j;
        tmp_2-=i;
        long long indexi=n1;
        indexi*=j;
        indexi=indexi+i;
        upper[tmp_2-1]=cormat[indexi];
//if((i==39001 &&j == 69999)||(i==1 && j==2))
 // printf("\n\n\n thread:  %f ",upper[tmp_2-1]," ",cormat[indexi]);
}

}


int CorMat_2(float* upper_tri, float * BOLD, int N, int L)
{
    long long M1 = (N-1); //computing the  size of correlaion matrix
    M1 *= N;
    M1 /= 2;
    long long total=N*N;//size of total correlation matrix

   // float * total_cormat = new float [total];

    preprocessing( BOLD,  N, L);//Preprocessing fMRI data in CPU
    
    cudaError_t cudaStat;
    cublasStatus_t stat;
    cublasHandle_t handle;
    stat = cublasCreate(&handle) ;
    if (stat != CUBLAS_STATUS_SUCCESS)
    {
        cout<<"Error in creating cublas handle";
        return stat;
    }
    
    
    float * devBOLD; //Allocating space in GPU for storing fMRI data
    cudaStat = cudaMalloc ((void**)&devBOLD, sizeof(float) * L * N) ;
    
    if (cudaStat != cudaSuccess)
    {
        cout<<"Error in Cuda Malloc";
        return cudaStat;
    }
    
    stat = cublasSetMatrix(N, L, sizeof(float), BOLD, N, devBOLD, N);//Copying fMRI data from CPU to GPU
    if (stat != CUBLAS_STATUS_SUCCESS)
     cout<<"Error in copying data to GPU";
    
    
    const float alpha = 1.0;
    const float beta = 0.0;
    
    float* devCormat;//allocating space in GPU for whole correlation matrix
    cudaMalloc ( (void**)&devCormat, sizeof(float) * total) ;
    
    stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, N,N,L,  &alpha, devBOLD, L, devBOLD, L, &beta, devCormat, N);//Performing matrix multiplication (fMRI data to its transpose)
    cudaDeviceSynchronize();
    
    if (stat != CUBLAS_STATUS_SUCCESS)
    {
        cout<<"Error performing multiplication";
        return stat;
    }
    
    float* dev_upper;//Allocating space for extracting upper triangle part
    cudaMalloc ( (void**)&dev_upper, sizeof(float) * M1) ;
    
    int block_size=1024;//number of threads
    long long grid_size=1+((total-1)/block_size);//number of blocks
    ker<<<grid_size,block_size>>>(devCormat,dev_upper,N,N);//performing kernel for extracting and reordering correlations from upper triangle
 memset((void*)upper_tri, 0, sizeof(float) *M1);
    gpuErrchk( cudaPeekAtLastError() );
    
    cudaMemcpy(upper_tri, dev_upper, sizeof(float) *M1, cudaMemcpyDeviceToHost);//copying upper triangle correlation matrix data back to CPU
    
    if (cudaStat != cudaSuccess)
    {
        cout<<"Error in Cuda memcpy back to cpu";
        return cudaStat;
    }
    
    cudaFree (devCormat);
    cudaFree (dev_upper);
    stat = cublasDestroy(handle);
    if (stat != CUBLAS_STATUS_SUCCESS)
    {
        cout<<"Error in destroy";
        return stat;
    }
    return 1;
    
    
}


////////////////////////////////////////

int CorMat_3(float* upper_tri, float * BOLD, int N, int L,long long OOO)
{
    //clock_t first,second;    
    size_t free;int ii=0;
    size_t total_mem;
    cudaMemGetInfo(&free,&total_mem);
    long long available_mem = free;
    available_mem/=sizeof(float);
    available_mem-=(N*L);//Getting available memory without
    
    int flag=1;

    preprocessing( BOLD,  N, L);//Preprocessing fMRI data

    cudaError_t cudaStat;
    cublasStatus_t stat;
    cublasHandle_t handle;
    long long upper_size=(N-1);//computinf size of total correlation matrix
    upper_size*=N;
    upper_size/=2;
    
    float * devBOLD;//initializing normalized fMRI data in GPU
    cudaStat = cudaMalloc ((void**)&devBOLD, sizeof(float) * L * N);

    stat = cublasSetMatrix(N, L, sizeof(float), BOLD, N, devBOLD, N);

    cublasCreate(&handle) ;
    
    const float alpha = 1.0;
    const float beta = 0.0;
    
    int block,N_prime;
    block=OOO;
    N_prime=N;
    
    float* add_uper_cpu=upper_tri;
    long long M1,temp,temp2=0,temp3=0;
    int so_far=0;
    int pak=0;
    float* devCormat;
    float* dev_upper;
    int ffl=0;
    long long old_cormat_fullsize;
    long long old_M1;
    long long cormat_fullsize;
    while(flag==1)
    {
        cout<<"this is block: "<<block<<"\n\n";        
        if(block==N_prime)//checking for the last chunk
           flag=0;
        
        temp = block;
        temp *= (block +1);
        temp /= 2;
        M1=N_prime;
        M1*=block;
        M1-=temp; //M1 is the size of upper triangle part of chunk


		if(pak!=0)
			{
			cudaFree (dev_upper);
			cudaFree (devCormat);

			}
            cormat_fullsize=block;
            cormat_fullsize*=N_prime;

            cudaStat=cudaMalloc ( (void**)&devCormat, sizeof(float) * cormat_fullsize) ;
            
            if (cudaStat != cudaSuccess)
            
            {
                cout<<"Error in Cuda Malloc and status is devcormat: "<<cudaStat;
                return cudaStat;
            }
            
            cudaStat =  cudaMalloc ( (void**)&dev_upper, sizeof(float) * M1) ;
            if (cudaStat != cudaSuccess)
            
            {
                cout<<"Error in Cuda Malloc and status is devcormat: "<<cudaStat;
                return cudaStat;
            }
  
            
            cout<<"\n IN PAK  0: "<<cormat_fullsize<<" " <<M1<<"*****";
            old_cormat_fullsize=cormat_fullsize;
            old_M1=M1;
            pak++;

        stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, block,N_prime,L,  &alpha, devBOLD+(so_far*L), L, devBOLD+(so_far*L), L, &beta, devCormat, block);//multiply block x L to L x N_prime = block x N_prime

        if (stat != CUBLAS_STATUS_SUCCESS)
        {
            cout<<"error in cublasSgemm, stat is: \n";
            cout<<stat<<"\n";
            return stat;
        }

        cudaDeviceSynchronize();

        temp2=block;
        temp2*=N_prime;

        int block_size=1024;
        long long grid_size=1+((temp2-1)/block_size);

        ker2<<<grid_size,block_size>>>(devCormat,dev_upper,block,N_prime,upper_size,N,ii,M1);
        
	memset((void*)add_uper_cpu, 0, sizeof(float) *M1); 
	
        cudaDeviceSynchronize();
        ii+=block;
        
        gpuErrchk( cudaPeekAtLastError() );
        
        cudaStat= cudaMemcpy(add_uper_cpu, dev_upper, sizeof(float) *M1, cudaMemcpyDeviceToHost);

        if (cudaStat != cudaSuccess)
        {
            cout<<"cudamalloc add_uper: \n";
            cout<<stat<<"\n";
            return stat;
        }

        temp3+=M1;
        add_uper_cpu=upper_tri+temp3;
        so_far+=block;
        
        if(N_prime>block)
        {
            N_prime=N_prime-block;
            block=remaining_N2( N_prime, L, available_mem);
          
            if(N_prime  <block)//checking last chunk
             block=N_prime;

        }
    }
    cudaFree (devBOLD);
    cudaFree (dev_upper);
    cudaFree (devCormat);

    stat = cublasDestroy(handle);
    if (stat != CUBLAS_STATUS_SUCCESS)
    {
        cout<<"error in destroy";
        return stat;
    }
    
    
    return 1;
}





///////////////////////////////////////


/*
 Copyright (C)Taban Eslami and Fahad Saeed
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

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
int CorMat_2(float* , float * , int , int );
int CorMat_3(float* , float * , int , int,long long);
long long remaining_N2(int N, int L,long long available_memory)
{
    long long x=available_memory;
    long long temp=N;
    temp*=2;
    x/=temp;
    return x;
}


long long remaining_N(int N, int L,int fflag)
{
        size_t free;
        size_t total_mem;
        cudaMemGetInfo(&free,&total_mem);
        long long available_mem = free;
        available_mem/=sizeof(float);
        long long NL=N;
        NL*=L;
        if (fflag==0)
        available_mem-=(NL);
        long long x=available_mem;
        long long temp=N;
        temp*=2;
        x/=temp;
        return x;
}


void preprocessing(float * BOLD_t, int N,int L)
{

      for (int i = 0; i < N; i++)
                {
                        float * row = BOLD_t + i * L;
                        double sum1 = 0, sum2 = 0;
                        for (int l = 0; l < L; l++)
                        {
                                sum1 += row[l];
                        }
                        sum1 /= L;
                        for (int l = 0; l < L; l++)
                        {
                                sum2 += (row[l] - sum1) * (row[l] - sum1);
                        }
                                   sum2 = sqrt(sum2);
                        for (int l = 0; l < L; l++)
                        {
                                if(sum2!=0)
                                row[l] = (row[l] - sum1) / sum2;
                                else
                                if(sum2==0)
                                row[l]=0;
                        }
                }

}



int main(int argc, char *argv[])
{
    int  k = 0, l = 0;
    int N,L;
    char wr;
    N = atoi(argv[1]);
    L = atoi(argv[2]);
    wr = *(argv[3]);
     cout<<"Number of voxels: "<<N<<"  "<<"Length of time series: "<<L<<"\n\n";
        srand(time(0));

    clock_t first,second;
    string name ="./";
    stringstream sstm;
    ifstream myReadFile;
    sstm.str("");
    sstm << name<<"data.txt";
    clock_t kho1,kho2;
    kho1=clock();
    string ad = sstm.str();
         myReadFile.open(ad.c_str());

        float * BOLD = new float [L * N];
        for (  k = 0; k < N; k++){
           for ( l = 0; l < L; l++)
        {
                myReadFile>>BOLD[k*L+l];//BOLD[l*N+k];

        }
        }

        myReadFile.close();
kho2=clock();
/////////////////////////////////////////////////////
        long long M11 = (N-1);
        M11 *= N;
        M11 /= 2;


cudaEvent_t start, stop;
float time;

float * upper_tri= new float [M11];

for(long long indii=0;indii<M11;indii++)
upper_tri[indii]=0;

long long OOO=remaining_N( N, L,0);

 long long temp_in=M11/2;
    long long temp_in2=M11;
     temp_in2*=3;
    temp_in2/=4;
    long long temp_in3=M11;
    temp_in3*=4;
    temp_in3/=5;


cout<<"\nComputing correlations ...";

if(N<=OOO)
{

cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord( start, 0 );

  //      cout<<"\n This means that we have enough memory ";

        first = clock();
        int u = CorMat_2(upper_tri,BOLD,N, L);

        second = clock();
        cout<<"\nRunning time for computing correlations: \n"<<(double)(second-first)/CLOCKS_PER_SEC<<" \n";
        delete []BOLD;
cudaEventRecord( stop, 0 );
cudaEventSynchronize( stop );

cudaEventElapsedTime( &time, start, stop );

}

if(N>OOO)
{

cudaEventCreate(&start);
cudaEventCreate(&stop);

cudaEventRecord( start, 0 );
first = clock();
CorMat_3(upper_tri,  BOLD,  N, L, OOO);
second = clock();
  cout<<"\nRunning time for computing correlations: \n"<<(double)(second-first)/CLOCKS_PER_SEC<<"\n ";


    delete []BOLD;
cudaEventRecord( stop, 0 );
cudaEventSynchronize( stop );

cudaEventElapsedTime( &time, start, stop );

}

if (wr == 'b'){
cout<<"\nWriting correlation values into the binary file ... \n";
ofstream OutFile;
OutFile.open("corrs.bin", ios::binary | ios::out);
OutFile.write((char*)upper_tri, sizeof(float)*M11);
OutFile.close();
cout<<"\nCorrelations are stored into the file corrs.bin \n";
}

if (wr == 't'){
cout<<"\nWriting correlation values into the text file ... \n";
ofstream correlations_print;
correlations_print.open("corrs.txt");
for(long long tab =0;tab<M11;tab++)

     {    
           correlations_print << upper_tri[tab] << '\n';
     }

correlations_print.close();
cout<<"\nCorrelations are stored into the text file corrs.txt \n";
}
return 0;

}

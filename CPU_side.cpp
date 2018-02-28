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
    N = atoi(argv[1]);
    L = atoi(argv[2]);
    cout<<N<<"  "<<L<<"\n\n";
        srand(time(0));

    clock_t first,second;
    string name ="/home/taban/Correlation_codes/code_wang/200000_500/";// "/home/taban/Correlation_codes/code_wang/100000_500/";
    stringstream sstm;
    ifstream myReadFile;
    sstm.str("");
    sstm << name<<"data1.txt";//"5row3column.txt";//"data1.txt";//"data_array_33001_215"<<".txt";//"data_array_18384_150"<<".txt";//"data_array_33001_215"<<".txt";//"data1"<<".txt";//"data_array_33001_215"<<".txt";
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
cout<<"\nfile reading: \n"<<(double)(kho2-kho1)/CLOCKS_PER_SEC<<"\n ======================= \n";
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




if(N<=OOO)
{

cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord( start, 0 );

        cout<<"\n This means that we have enough memory ";

        first = clock();
        int u = CorMat_2(upper_tri,BOLD,N, L);

        second = clock();
        cout<<"\ntime is: \n"<<(double)(second-first)/CLOCKS_PER_SEC<<"\n ======================= \n";
        delete []BOLD;
cudaEventRecord( stop, 0 );
cudaEventSynchronize( stop );

cudaEventElapsedTime( &time, start, stop );

}

if(N>OOO)
{
cout<<"\n this is OOO: "<<OOO<<"\n";

cudaEventCreate(&start);
cudaEventCreate(&stop);

cudaEventRecord( start, 0 );
first = clock();
CorMat_3(upper_tri,  BOLD,  N, L, OOO);
second = clock();
  cout<<"\ntime is: \n"<<(double)(second-first)/CLOCKS_PER_SEC<<"\n ======================= \n";


    delete []BOLD;
cudaEventRecord( stop, 0 );
cudaEventSynchronize( stop );

cudaEventElapsedTime( &time, start, stop );

}



cout<<"\n event time is : "<<time/1000<<"\n";

//cout<<"\n data: - \n"<<temp_in<<": "<<upper_tri[temp_in]<<" -  "<<temp_in2<<": "<<upper_tri[temp_in2]<<" -  "<<temp_in3<<": "<<upper_tri[temp_in3]<<"\n "<<upper_tri[0]<<"  "<<upper_tri[1]<<"  "<<upper_tri[100]<<"  "<<upper_tri[1233]<<"\n "<<upper_tri[M11-100]<<" "<<upper_tri[M11-1]<<" "<<upper_tri[M11-2]<<"\n";
cout<<"\n data: - \n"<<upper_tri[temp_in]<<" "<<upper_tri[temp_in2]<<" "<<upper_tri[temp_in3]<<"\n "<<upper_tri[0]<<" "<<upper_tri[1]<<" "<<upper_tri[100]<<" "<<upper_tri[1233]<<"\n "<<upper_tri[M11-100]<<" "<<upper_tri[M11-1]<<" "<<upper_tri[M11-2]<<"\n";
//cudaEventDestroy( start );
//cudaEventDestroy( stop );
      cout<<"\n --------------------------- \n";

cout<<upper_tri[656767]<<"\n\n";

return 0;

}

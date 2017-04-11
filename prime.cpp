#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

using namespace std;

int CountPrimes(unsigned long long int n0, unsigned long long int n1, unsigned long int *basePrimes, int nBasePrimes)
{
	// Use the sieve of Eratosthenes to count the number of primes between n0 and n1-1 inclusive
	// Assumptions:
	// * basePrimes is a pointer to a list of known primes in ascending order
	// n0 > max(basePrimes)
	// * n1 <= max(basePrimes)^2 (a composite number must have a prime factor <= its square root)
	// * n1-n0 is small enough so we can allocate an array of (n1-n0)*64 bits.
	int nPrimes = 0;
	int i;
	unsigned long long *working;
	working = (unsigned long long int*)malloc((n1 - n0) * sizeof(unsigned long long int));
	memset(working, 0, (n1 - n0) * sizeof(unsigned long long int));
	for (i = 0; i < nBasePrimes; i++)
	{
		// find starting point
		long long int i0;
		if ((n0 % basePrimes[i]) == 0)
			i0 = n0;
		else
			i0 = (n0 / basePrimes[i] + 1)*basePrimes[i];
		while (i0 < n1)
		{
			working[i0 - n0] = 1;
			i0 += basePrimes[i];
		}
	}
	for (i = 0; i < (n1 - n0); i++)
		if (!working[i])
			nPrimes++;
	free(working);
	return nPrimes;
}           
int main(int argc, char** argv)
{
	int rank, comm_size;
	unsigned long long int sum=0;
	double elapsed_time;
	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time=-MPI_Wtime();
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	
	unsigned long long int nPrimes;
	unsigned long long int nPrimes_divide;
  int i, j;
    
#define MAX_PRIME 100000 // we are counting primes up to 1e10, so we need a base primes list up to sqrt(1e10) = 1e5
#define SQRT_MAX_PRIME 317 // the largest number we will need to cross out multiples of for our base prime list

	unsigned long int *working = (unsigned long int *)malloc((MAX_PRIME + 1) * sizeof(unsigned long int));
	unsigned long int *primes = (unsigned long int *)malloc(MAX_PRIME * sizeof(unsigned long int));
	unsigned long int *primes_divide = (unsigned long int *)malloc(MAX_PRIME * sizeof(unsigned long int));

		
   

	// construct the base prime list (primes up to 1e5) using the sieve of Eratosthenes
	

for ( i = 0; i < MAX_PRIME; i++ )
	{
		working[i] = 1;
	}
	working[0] = working[1] = 0;

	// sieve
for ( i = 2; i < SQRT_MAX_PRIME; i++ )
	{
		if ( working[i] )
		{
			for ( j = 2*i; j < MAX_PRIME; j += i )
			{
				working[j] = 0;
			}
		}
	}
    //saving the base prime list twice so it can be divided later from the sum as it was calculated on each rank
    nPrimes = 0;
for ( i = 0; i < MAX_PRIME; i++ )
		if ( working[i] )
			primes[nPrimes++] = i; 

    nPrimes_divide = 0;
for ( i = 0; i < MAX_PRIME; i++ )
		if ( working[i] )
			primes_divide[nPrimes_divide++] = i;

	unsigned long long int ip,  ip_data;
	unsigned long long int delta = 100000000;
	unsigned long long int pmax = 10000000000;
	int nBase = nPrimes;
  ip = primes[nPrimes-1] + 1; 
	ip_data=ip;
  unsigned long long int size = (pmax / delta) +1 ;
	unsigned long long int localdata_size = size/comm_size; 
	unsigned long long int *data;
	unsigned long long int *data_rem;
	unsigned long long int *localdata;
	unsigned long long int last_data;
	unsigned long long int *nPrimes_sum;
	int divided = size/comm_size;
	int new_size = size; 
    long long int rem_size = new_size-((divided*comm_size)-1);
	data_rem = (unsigned long long int*)malloc(sizeof(unsigned long long int) * size);//remaining data of the data array after divided by the comm_size
	data = (unsigned long long int*)malloc(sizeof(unsigned long long int) * size); //declare array of data
	localdata = (unsigned long long int*)malloc(sizeof(unsigned long long int) * size); //declare array of data
	nPrimes_sum = (unsigned long long int*)malloc(sizeof(unsigned long long int) * size);//sum of the nPrimes
	data[0]=ip; //first element of data array

if(rank==0)
    { 
    //populating the data array on rank 0
while ( ip_data < pmax )
	{
for (int i = 1; i < size; i++) {
			unsigned long long int next = ip_data + delta;
		if (next > pmax)
			next = pmax;
		ip_data = next;		
			data[i] = ip_data;						
			}
		}
    }
	last_data= ((divided*comm_size-1)*delta)+ data[0]; //last data if the data array was equally divided between the ranks
	data_rem[0]=last_data; //first element of the data_rem array
   
if (rank==0)
    {
        printf("data_rem[0]: %lld\n", data_rem[0]);	
	//populating data_rem array
	for (int i = 1; i < rem_size; i++) {		

			unsigned long long int next = last_data + delta;
		if (next > pmax)
			next = pmax;
		last_data = next;
		
			data_rem[i] = last_data;						
			} 	
		printf("data_rem_last: %lld\n", data_rem[rem_size-1]);
    }

MPI_Scatter(data, localdata_size, MPI_LONG_LONG_INT, localdata, localdata_size, MPI_LONG_LONG_INT, 0, 
	MPI_COMM_WORLD); //scattering the data array to every rank storing it in the localdata array

	unsigned long long int next_local ; 	
for (i = 0; i < size; i++){ 
		if (next_local > pmax)
			next_local = pmax;
	ip = localdata[i];
if (localdata[i] != localdata[divided]){
    next_local = ip + delta;
			if (next_local > pmax)
			next_local = pmax;
			printf("counting on processor %d from %lld to %lld\n", rank, ip, next_local);
			nPrimes += CountPrimes(ip, next_local, primes, nBase);
			ip = next_local;
	}
    }
if (rank ==0){
	unsigned long long int ip_rem = data_rem[1];
if(data_rem[0] < pmax)
	{
	
while (ip_rem < pmax)
	{
		unsigned long long int next_rem = ip_rem + delta;
		if (next_rem > pmax)
			next_rem = pmax;
		printf("counting on %d from %lld to %lld\n",rank, ip_rem, next_rem);
		nPrimes += CountPrimes(ip_rem, next_rem, primes, nBase);
		ip_rem = next_rem;	
	}    
    }
    }
if (rank != 0)
    { 
    nPrimes = nPrimes - nPrimes_divide;
   }
    printf("prime counting on %d function(%lld) = %lld\n",rank, pmax, nPrimes);
    MPI_Gather(&nPrimes, 1, MPI_LONG_LONG_INT, nPrimes_sum,1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD); //gathering the number of primes from each rank
    MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime(); //stopping the timer
if ( rank ==0){
    //calculating the sum of the primes
    for (int j = 0; j< comm_size; j++){
    sum = sum + nPrimes_sum[j];
    }
    printf("prime counting function(%lld)  = %lld\n", pmax, sum);
    printf("Total elapsed time: %10.6fs\n", elapsed_time);
    }
    free(primes);
	MPI_Finalize();
	return 0;
    }

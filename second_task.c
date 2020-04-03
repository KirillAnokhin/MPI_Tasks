#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define ROOT 0
#define BCAST 1
#define REDUCE 2
#define SCATTER 3
#define GATHER 4

void timer(int rank, double *arr) {
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == ROOT)
		arr[0] = MPI_Wtime();
}

void run_function(int n_func) {
	int bcast = 0;
	int sendbuf, recvbuf;
	int size = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int *sbuf = calloc(size, sizeof(int));
	int *rbuf = calloc(size, sizeof(int));

	if(n_func == BCAST)
		MPI_Bcast(&bcast, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	else if(n_func == REDUCE)
		MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_PROD, ROOT, MPI_COMM_WORLD);
	else if(n_func == SCATTER)
		MPI_Scatter(sbuf, 1, MPI_INT, rbuf, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	else if(n_func == GATHER)
		MPI_Gather(sbuf, 1, MPI_INT, rbuf, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	free(sbuf);
	free(rbuf);
}

double calculate_average (int rank, int n_func) {
	double sigma = 1.;
	int n_runs = 50;
	double avrg = 0;
	double precision = MPI_Wtick();
	precision *= 100;
	if (rank == ROOT) {
		while(sigma > precision) {
			MPI_Bcast(&n_runs, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
			double *starts = calloc(n_runs, sizeof(double));
			double *ends = calloc(n_runs, sizeof(double));
			for(int i = 0; i < n_runs; i++) {
				timer(rank, starts + i);
				run_function(n_func);
				timer(rank, ends + i);
				ends[i] -= starts[i];
			}

			for(int i = 0; i < n_runs; i++) 
				avrg += ends[i];

			avrg /= n_runs;
			double dispertion  = 0;
			for(int i = 0; i < n_runs; i++)
				dispertion += ((ends[i] - avrg) * (ends[i] - avrg));
			sigma = sqrt(dispertion/n_runs);
			printf("dispertion %0.9lf sigma: %0.9lf precision: %0.9lf n_runs: %d \n", 
				       dispertion, sigma, precision, n_runs);
			free(starts);
			free(ends);
			n_runs *= 2;
		}
		n_runs = 0;
		MPI_Bcast(&n_runs, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	} else {
		while(n_runs > 0) {
			MPI_Bcast(&n_runs, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
			for(int i = 0; i < n_runs; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				run_function(n_func);
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}

	}
	return avrg;

}

int main (int argc, char *argv[]) 
{
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double bcast_time = calculate_average(rank, BCAST);
	double reduce_time = calculate_average(rank, REDUCE);
	double scatter_time = calculate_average(rank, SCATTER);
	double gather_time = calculate_average(rank, GATHER);
	if(rank == ROOT) {
		printf("Main:\n");
		printf("bcast_time: %0.9lf\nreduce_time: %0.9lf\nscatter_time: %0.9lf\ngather_time: %0.9lf\n", 
				bcast_time, reduce_time, scatter_time, gather_time);
	}
	MPI_Finalize();
	return 0;
}

#include <stdio.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char *argv[])
{
	int ProcNum, ProcRank, RecvRank;

	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	if(ProcRank == 0) {
		printf("Hello from proc %d\n", ProcRank);
		for(int i = 1; i < ProcNum; i++) {
			MPI_Send(&ProcRank, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Recv(&RecvRank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &Status);
		}
	}
	else {
		MPI_Recv(&RecvRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		printf("Hello from proc %d\n", ProcRank);
		MPI_Send(&ProcRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}

#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	int ProcNum, ProcRank, RecvRank;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	if(ProcRank != 0) 
		MPI_Recv(&RecvRank, 1, MPI_INT, ProcRank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status); 

	printf("Hello from process %d\n", ProcRank);

	if(ProcRank != ProcNum - 1)
		MPI_Send(&ProcRank, 1, MPI_INT, ProcRank + 1, 0, MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

#include <macros.hpp>

#ifdef COSMO_MPI
#include <mpi.h>

bool IS_PARALLEL()
{
    int f;
    MPI_Initialized(&f);
    return f;
}

int CURRENT_PROCESS()
{
    if(!IS_PARALLEL())
        return 0;

    int p;
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    return p;
}

int NUM_PROCESSES()
{
    if(!IS_PARALLEL())
        return 1;

    int n;
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    return n;
}

#else
bool IS_PARALLEL() { return false; }
int CURRENT_PROCESS() { return 0; }
int NUM_PROCESSES() { return 1; }
#endif

#ifdef COSMO_OMP
#include <omp.h>

int TOTAL_NUM_THREADS() { return omp_get_num_threads(); }
int CURRENT_THREAD_NUM() { return omp_get_thread_num(); }
#else
int TOTAL_NUM_THREADS() { return 1; }
int CURRENT_THREAD_NUM() { return 0; }
#endif


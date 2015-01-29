#ifndef COSMO_PP_COSMO_MPI_HPP
#define COSMO_PP_COSMO_MPI_HPP

#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <macros.hpp>

class CosmoMPI
{
private:
    CosmoMPI()
    {
#ifdef COSMO_MPI
#ifdef CHECKS_ON
        int hasMpiInitialized;
        MPI_Initialized(&hasMpiInitialized);
        check(!hasMpiInitialized, "MPI already initialized");
#endif

        MPI_Init(NULL, NULL);
#endif
    }

    ~CosmoMPI()
    {
#ifdef COSMO_MPI
#ifdef CHECKS_ON
        int hasMpiInitialized;
        MPI_Initialized(&hasMpiInitialized);
        check(hasMpiInitialized, "MPI not initialized");

        int hasMpiFinalized;
        MPI_Finalized(&hasMpiFinalized);
        check(!hasMpiFinalized, "MPI already finalized");
#endif
        MPI_Finalize();
#endif
    }
public:
	static CosmoMPI& create()
    {
        static CosmoMPI c;
        return c;
    }

    int processId() const
    {
#ifdef COSMO_MPI
#ifdef CHECKS_ON
        int hasMpiInitialized;
        MPI_Initialized(&hasMpiInitialized);
        check(hasMpiInitialized, "MPI not initialized");

        int hasMpiFinalized;
        MPI_Finalized(&hasMpiFinalized);
        check(!hasMpiFinalized, "MPI already finalized");
#endif
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
#else
        return 0;
#endif
    }

    int numProcesses() const
    {
#ifdef COSMO_MPI
#ifdef CHECKS_ON
        int hasMpiInitialized;
        MPI_Initialized(&hasMpiInitialized);
        check(hasMpiInitialized, "MPI not initialized");

        int hasMpiFinalized;
        MPI_Finalized(&hasMpiFinalized);
        check(!hasMpiFinalized, "MPI already finalized");
#endif
        int n;
        MPI_Comm_size(MPI_COMM_WORLD, &n);
        return n;
#else
        return 1;
#endif
    }

    bool isMaster() const
    {
        return (processId() == 0);
    }
};

#endif


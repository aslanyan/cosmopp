#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <cosmo_mpi.hpp>
#include <macros.hpp>

CosmoMPI::CosmoMPI()
{
#ifdef COSMO_MPI
#ifdef CHECKS_ON
    int hasMpiInitialized;
    MPI_Initialized(&hasMpiInitialized);
    check(!hasMpiInitialized, "MPI already initialized");
#endif

    MPI_Init(NULL, NULL);
#endif
    commTag_ = 1000;
}

CosmoMPI::~CosmoMPI()
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

int
CosmoMPI::processId() const
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

int
CosmoMPI::numProcesses() const
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

void
CosmoMPI::barrier() const
{
#ifdef COSMO_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

int
CosmoMPI::getCommTag()
{
    barrier();

    commTag_ += 10 * numProcesses();

    return commTag_;
}

int
CosmoMPI::send(int dest, void *buf, int count, DataType type, int tag)
{
    check(dest >= 0 && dest < numProcesses(), "invalid destination " << dest);
    check(dest != processId(), "cannot send to myself");

    check(type >= 0 && type < DATA_TYPE_MAX, "");

#ifdef COSMO_MPI
    MPI_Datatype mpiType;
    switch(type)
    {
    case DOUBLE:
        mpiType = MPI_DOUBLE;
        break;
    case INT:
        mpiType = MPI_INT;
        break;
    case LONG:
        mpiType = MPI_LONG;
        break;
    default:
        check(false, "");
        break;
    }

    const int res = MPI_Send(buf, count, mpiType, dest, tag, MPI_COMM_WORLD);
    check(res == MPI_SUCCESS, "send failed");
    return res;
#else
    check(false, "");
    return 0;
#endif
}

int
CosmoMPI::recv(int source, void *buf, int count, DataType type, int tag)
{
    check(source >= 0 && source < numProcesses(), "invalid source " << source);
    check(source != processId(), "cannot receive from myself");

    check(type >= 0 && type < DATA_TYPE_MAX, "");

#ifdef COSMO_MPI
    MPI_Datatype mpiType;
    switch(type)
    {
    case DOUBLE:
        mpiType = MPI_DOUBLE;
        break;
    case INT:
        mpiType = MPI_INT;
        break;
    case LONG:
        mpiType = MPI_LONG;
        break;
    default:
        check(false, "");
        break;
    }

    MPI_Status s;
    const int res = MPI_Recv(buf, count, mpiType, source, tag, MPI_COMM_WORLD, &s);
    check(res == MPI_SUCCESS, "receive failed");
    return res;
#else
    check(false, "");
    return 0;
#endif
}


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

#ifdef COSMO_MPI
namespace
{

MPI_Datatype typeToMPIType(CosmoMPI::DataType type)
{
    MPI_Datatype mpiType;
    switch(type)
    {
    case CosmoMPI::DOUBLE:
        mpiType = MPI_DOUBLE;
        break;
    case CosmoMPI::INT:
        mpiType = MPI_INT;
        break;
    case CosmoMPI::LONG:
        mpiType = MPI_LONG;
        break;
    default:
        check(false, "");
        break;
    }
    return mpiType;
}

MPI_Op opToMPIOp(CosmoMPI::ReduceOp op)
{
    MPI_Op mpiOp;
    switch(op)
    {
        case CosmoMPI::SUM:
            mpiOp = MPI_SUM;
            break;
        case CosmoMPI::MAX:
            mpiOp = MPI_MAX;
            break;
        case CosmoMPI::MIN:
            mpiOp = MPI_MIN;
            break;
        case CosmoMPI::PROD:
            mpiOp = MPI_PROD;
            break;
        default:
            check(false, "");
            break;
    }
    return mpiOp;
}

} // namespace
#endif

int
CosmoMPI::send(int dest, void *buf, int count, DataType type, int tag)
{
    check(dest >= 0 && dest < numProcesses(), "invalid destination " << dest);
    check(dest != processId(), "cannot send to myself");

    check(type >= 0 && type < DATA_TYPE_MAX, "");

#ifdef COSMO_MPI
    const MPI_Datatype mpiType = typeToMPIType(type);
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
    const MPI_Datatype mpiType = typeToMPIType(type);
    MPI_Status s;
    const int res = MPI_Recv(buf, count, mpiType, source, tag, MPI_COMM_WORLD, &s);
    check(res == MPI_SUCCESS, "receive failed");
    return res;
#else
    check(false, "");
    return 0;
#endif
}

int
CosmoMPI::reduce(void *send, void *recv, int count, DataType type, ReduceOp op)
{
    check(type >= 0 && type < DATA_TYPE_MAX, "");
    check(op >= 0 && op < REDUCE_OP_MAX, "");

#ifdef COSMO_MPI
    const MPI_Datatype mpiType = typeToMPIType(type);
    const MPI_Op mpiOp = opToMPIOp(op);
    const int res = MPI_Reduce(send, recv, count, mpiType, mpiOp, 0, MPI_COMM_WORLD);
    check(res == MPI_SUCCESS, "reduce failed");
    return res;
#else
    check(false, "");
    return 0;
#endif
}

int
CosmoMPI::bcast(void *data, int count, DataType type)
{
    check(type >= 0 && type < DATA_TYPE_MAX, "");

#ifdef COSMO_MPI
    const MPI_Datatype mpiType = typeToMPIType(type);
    const int res = MPI_Bcast(data, count, mpiType, 0, MPI_COMM_WORLD);
    check(res == MPI_SUCCESS, "bcast failed");
    return res;
#else
    return 0;
#endif
}


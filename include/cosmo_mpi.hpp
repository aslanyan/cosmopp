#ifndef COSMO_PP_COSMO_MPI_HPP
#define COSMO_PP_COSMO_MPI_HPP

class CosmoMPI
{
private:
    CosmoMPI();
    ~CosmoMPI();

public:
	static CosmoMPI& create()
    {
        static CosmoMPI c;
        return c;
    }

    int processId() const;
    int numProcesses() const;
    bool isMaster() const { return (processId() == 0); }
    void barrier() const;
    int getCommTag();

    enum DataType { DOUBLE = 0, INT, LONG, DATA_TYPE_MAX };
    enum ReduceOp { SUM = 0, MAX, MIN, PROD, REDUCE_OP_MAX };

    int send(int dest, void *buf, int count, DataType type, int tag);
    int recv(int source, void *buf, int count, DataType type, int tag);

    int reduce(void *send, void *recv, int count, DataType type, ReduceOp op);

private:
    int commTag_;
};

#endif


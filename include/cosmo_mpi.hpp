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

private:
    int commTag_;
};

#endif


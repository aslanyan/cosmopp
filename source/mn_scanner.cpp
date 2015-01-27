#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <mn_scanner.hpp>

#include <multinest.h>


MnScanner::MnScanner(int nPar, Math::LikelihoodFunction& like, int nLive, std::string fileRoot, bool accurateEvidence) : n_(nPar), like_(like), nLive_(nLive), paramsStarting_(nPar, 0), paramNames_(nPar), paramsBest_(nPar, 0), paramsMean_(nPar, 0), paramsStd_(nPar, 0), paramsCurrent_(nPar, 0), paramPriors_(nPar), paramsFixed_(nPar, 0), fileRoot_(fileRoot), accurateEvidence_(accurateEvidence)
{
}

void
MnScanner::setParam(int i, const std::string& name, double min, double max)
{
    check(i >= 0 && i < n_, "invalid index " << i);
    check(max >= min, "");

    if(min == max)
    {
        setParamFixed(i, name, min);
        return;
    }

    paramNames_[i] = name;
    paramPriors_[i].clear();

    const double epsilon = 1e-7; // for numerics
    paramPriors_[i][0 - epsilon] = min;
    paramPriors_[i][1 + epsilon] = max;
}

void
MnScanner::setParamFixed(int i, const std::string& name, double val)
{
    check(i >= 0 && i < n_, "invalid index " << i);
    paramNames_[i] = name;
    paramPriors_[i].clear();
    paramsFixed_[i] = val;
}

void
MnScanner::setParamGauss(int i, const std::string& name, double mean, double sigma)
{
    check(i >= 0 && i < n_, "invalid index " << i);
    check(sigma >= 0, "");

    if(sigma == 0)
    {
        setParamFixed(i, name, mean);
        return;
    }

    paramNames_[i] = name;

    const double xMin = mean - 10 * sigma, xMax = mean + 10 * sigma;
    const int N = 100000;

    const double epsilon = 1e-7; // for numerics
    paramPriors_[i][0 - epsilon] = xMin;
    double xPrev = xMin;
    double yPrev = 0;
    double res = 0;
    const double delta = (xMax - xMin) / N;

    const double norm = 1.0 / (std::sqrt(2 * Math::pi) * sigma);
    for(int j = 1; j <= N; ++j)
    {
        const double x = xMin + j * delta;
        const double y = norm * std::exp(-(x - mean) * (x - mean) / (2 * sigma * sigma));

        res += (x - xPrev) * (y + yPrev) / 2.0;
        
        check(res <= 1.0, "");

        paramPriors_[i][res] = x;
        xPrev = x;
        yPrev = y;
    }

    paramPriors_[i][1 + epsilon] = xMax;
}

void
MnScanner::setParamGeneral(int i, const std::string& name, double min, double max, const Math::RealFunction& distrib)
{
    check(i >= 0 && i < n_, "invalid index " << i);
    check(max > min, "");

    paramNames_[i] = name;
    paramPriors_[i].clear();
    const int N = 100000;

    const double epsilon = 1e-7; // for numerics
    paramPriors_[i][0 - epsilon] = min;
    double xPrev = min;
    double yPrev = 0;
    double res = 0;
    const double delta = (max - min) / N;

    for(int j = 1; j <= N; ++j)
    {
        const double x = min + j * delta;
        const double y = distrib.evaluate(x);

        res += (x - xPrev) * (y + yPrev) / 2.0;
        
        xPrev = x;
        yPrev = y;
    }

    const double norm = 1.0 / res;
    xPrev = min;
    yPrev = 0;
    res = 0;
    for(int j = 1; j <= N; ++j)
    {
        const double x = min + j * delta;
        const double y = norm * distrib.evaluate(x);

        res += (x - xPrev) * (y + yPrev) / 2.0;
        check(res <= 1.0 + 1e-5, "");

        paramPriors_[i][res] = x;
        xPrev = x;
        yPrev = y;
    }

    paramPriors_[i][1 + epsilon] = max;
}

void myLogLike(double *cube, int &ndim, int &npars, double &lnew, void *context)
{
    MnScanner* scanner = (MnScanner*) context;
    scanner->logLike(cube, ndim, npars, lnew);
}

void
MnScanner::logLike(double *cube, int &ndim, int &npars, double &lnew)
{
    check(ndim == n_ - nFixed_, "");
    check(npars == n_ - nFixed_, "");

    int j = 0;

    for(int i = 0; i < n_; ++i)
    {
        double x;
        if(paramPriors_[i].empty())
        {
            x = paramsFixed_[i];
        }
        else
        {
            x = paramPriors_[i].evaluate(cube[j]);
            cube[j] = x;
            ++j;
        }
        paramsCurrent_[i] = x;
    }

    lnew = - like_.calculate(&(paramsCurrent_[0]), n_) / 2.0;
}

void myDumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSLogZ, double &logZerr, void *context)
{
    MnScanner* scanner = (MnScanner*) context;
    scanner->dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, INSLogZ, logZerr);
}

void
MnScanner::dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSLogZ, double &logZerr)
{
    StandardException exc;
    check(nPar == n_ - nFixed_, "");

	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;


    std::stringstream postFileName;
    postFileName << fileRoot_ << "posterior.txt";
    std::ofstream outPost(postFileName.str().c_str());

    if(!outPost)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into posterior distribution file " << postFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];

    for(i = 0; i < nSamples; ++i)
    {
        for(j = 0; j < nPar + 1; ++j)
            outPost << postdist[i][j] << ' ';
        outPost << postdist[i][nPar + 1] << std::endl;
    }
    outPost.close();
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];


    std::stringstream constrFileName;
    constrFileName << fileRoot_ << "parameter_constraints.txt";
    std::ofstream outConstr(constrFileName.str().c_str());
    
    if(!outConstr)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into constraints file " << constrFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }


    outConstr << "Constraints from " << nSamples << " samples:" << std::endl;
    outConstr << "log(Z) = " << logZ << "+-" << logZerr << std::endl;
    outConstr << "INS log(Z) = " << INSLogZ << "+-" << logZerr << std::endl;
	// parameter constraints
    j = 0;
    for(i = 0; i < n_; ++i)
    {
        if(paramPriors_[i].empty())
        {
            paramsMean_[i] = paramsFixed_[i];
            paramsStd_[i] = 0.0;
            paramsBest_[i] = paramsFixed_[i];
        }
        else
        {
            paramsMean_[i] = paramConstr[0][j];
            paramsStd_[i] = paramConstr[0][nPar + j];
            paramsBest_[i] = paramConstr[0][2 * nPar + j];
            ++j;
        }
        
        outConstr << paramNames_[i] << ":   " << paramsBest_[i] << "    " << paramsMean_[i] << "+-" << paramsStd_[i] << std::endl;
    }
    outConstr.close();
}

void
MnScanner::run(bool res)
{
    StandardException exc;

    nFixed_ = 0;
    for(int i = 0; i < n_; ++i)
    {
        if(paramPriors_[i].empty())
            ++nFixed_;
    }

    check(nFixed_ < n_, "cannot have all of the parameters fixed");

	int IS = 1;					// do Nested Importance Sampling?
	int mmodal = 0;					// do mode separation?
	int ceff = (accurateEvidence_ ? 0 : 1); 	// run in constant efficiency mode?
	double efr = (accurateEvidence_ ? 0.3 : 0.8);  // set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = n_ - nFixed_;					// dimensionality (no. of free parameters)
	int nPar = n_ - nFixed_;					// total no. of parameters including free & derived parameters
	int nClsPar = n_ - nFixed_;				// no. of parameters to do mode separation on
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	//char root[100] = fileRoot_.c_str();			// root for output files
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					// need feedback on standard output?
	int resume = res;					// resume from a previous job?
	int outfile = 1;				// write output files?

	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
#ifdef COSMO_MPI
    int hasMpiInit;
    MPI_Initialized(&hasMpiInit);
    if(hasMpiInit)
        initMPI = 0;
#endif

	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
    void* context = (void*) this;

    // Creating the paramnames file
    std::stringstream paramNamesFileName;
    paramNamesFileName << fileRoot_ << ".paramnames";
    std::ofstream outPar(paramNamesFileName.str().c_str());

    if(!outPar)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into paramnames file " << paramNamesFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    for(int i = 0; i < n_; ++i)
    {
        outPar << paramNames_[i] << '\t' << paramNames_[i] << std::endl;
    }
    outPar.close();

	// calling MultiNest

    try
    {
	    nested::run(IS, mmodal, ceff, nLive_, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, fileRoot_.c_str(), seed, pWrap, fb, resume, outfile, initMPI,
	        logZero, maxiter, myLogLike, myDumper, context);
    } 
    catch (std::exception& e)
    {
        dumpInfo(e.what());
        throw e;
    }
}

void
MnScanner::dumpInfo(const char* error)
{
    std::ofstream out("mn_scanner_error_log.txt");
    out << "Multinest Scanner ran into a problem:" << std::endl << error << std::endl << "The values of the parameters:" << std::endl;
    for(int i = 0; i < n_; ++i)
    {
        out << paramNames_[i] << " = " << std::setprecision(30) << paramsCurrent_[i] << std::endl;
    }
    out.close();
}

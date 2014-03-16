#ifndef COSMO_PP_BEST_FIT_HPP
#define COSMO_PP_BEST_FIT_HPP

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <likelihood_function.hpp>

#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"

namespace Math
{

class BestFitMinuitCalculator : public ROOT::Minuit2::FCNBase
{
public:
    BestFitMinuitCalculator(int nPar, LikelihoodFunction& like, std::string fileRoot, unsigned long* calls = NULL) : nPar_(nPar), like_(like), calls_(calls), fileRoot_(fileRoot)
    {
        check(nPar > 0, "Invalid number of parameters " << nPar << ". Must be positive.");

        std::stringstream fileName;
        fileName << fileRoot_ << ".txt";
        std::ofstream out(fileName.str().c_str());
        StandardException exc;
        if(!out)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into file " << fileName.str() << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }
        out.close();
    }

    ~BestFitMinuitCalculator()
    {
    }

    virtual double operator()(const std::vector<double>& par) const
    {
        check(par.size() == nPar_, "the number of parameters must be " << nPar_ << ", however " << par.size() << " provided");
        std::vector<double> parCopy(par);

        std::stringstream fileName;
        fileName << fileRoot_ << ".txt";
        std::ofstream out(fileName.str().c_str(), std::ios::app);
        StandardException exc;
        if(!out)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into file " << fileName.str() << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        const double res = like_.calculate(&(parCopy[0]), nPar_);
        out << "1 " << res;
        for(int i = 0; i < nPar_; ++i)
            out << ' ' << par[i];
        out << std::endl;
        out.close();

        if(calls_)
            ++(*calls_);

        return res;
    }

    unsigned long calls() const { return (calls_ ? *calls_ : 0); }

    double Up() const {return 1.;}

private:
    int nPar_;
    Math::LikelihoodFunction& like_;
    unsigned long *calls_;
    std::string fileRoot_;
};

class BestFit
{
public:
    BestFit(int nPar, LikelihoodFunction& like, std::string fileRoot) : calculator_(nPar, like, fileRoot, &calls_), calls_(0), nPar_(nPar), fileRoot_(fileRoot), paramNames_(nPar), starting_(nPar, 0.0), min_(nPar, 0.0), max_(nPar, 1.0), error_(nPar, 0.01)
    {
        check(nPar > 0, "Invalid number of parameters " << nPar << ". Must be positive.");

        for(int i = 0; i < nPar; ++i)
        {
            std::stringstream name;
            name << "Parameter_" << i;
            paramNames_[i] = name.str();
        }
    }

    void setParam(int i, const std::string& name, double min, double max, double starting = std::numeric_limits<double>::max(), double accuracy = 0.0)
    {
        check(i >= 0 && i < nPar_, "invalid index " << i);
        check(max > min, "max = " << max << ", min = " << min << ". need to have max > min.");
        paramNames_[i] = name;
        min_[i] = min;
        max_[i] = max;
        starting_[i] = (starting == std::numeric_limits<double>::max() ? (max + min) / 2.0 : starting);
        check(starting_[i] >= min && starting_[i] <= max, "invalid starting value " << starting);
        error_[i] = (accuracy == 0.0 ? (max - min) / 100 : accuracy);
        check(error_[i] > 0, "invalid accuracy " << accuracy);
    }
    
    void setParamGauss(int i, const std::string& name, double mean, double sigma, double starting = std::numeric_limits<double>::max(), double accuracy = 0.0)
    {
        setParam(i, name, mean - 5 * sigma, mean + 5 * sigma, starting, accuracy);
    }

    void run()
    {
        ROOT::Minuit2::MnUserParameters upar;
        for(int i = 0; i < nPar_; ++i)
            upar.Add(paramNames_[i].c_str(), starting_[i], error_[i], min_[i], max_[i]);

        ROOT::Minuit2::MnMigrad migrad(calculator_, upar);
        ROOT::Minuit2::FunctionMinimum minRes = migrad();
        ROOT::Minuit2::MnUserParameters result = minRes.UserParameters();

        std::stringstream fileName;
        fileName << fileRoot_ << "best_fit.txt";
        std::ofstream out(fileName.str().c_str());
        if(!out)
        {
            StandardException exc;
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into output file " << fileName.str() << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }
        out << "Total number of calls: " << calculator_.calls() << std::endl;
        std::vector<double> params(nPar_);
        for(int i = 0; i < nPar_; ++i)
            params[i] = result.Value(paramNames_[i].c_str());

        const double bestFitLike = calculator_(params);
        out << "Best fit likelihood = " << bestFitLike << std::endl;
        for(int i = 0; i < nPar_; ++i)
            out << paramNames_[i] << ":\t" << params[i] << " +/- " << result.Error(paramNames_[i].c_str()) << std::endl;
        out.close();
    }

private:
    unsigned long calls_;
    BestFitMinuitCalculator calculator_;
    const int nPar_;
    const std::string fileRoot_;
    std::vector<std::string> paramNames_;
    std::vector<double> starting_, min_, max_, error_;
};

} // namespace Math

#endif


#ifndef COSMO_CPP_SCALE_FACTOR_HPP
#define COSMO_CPP_SCALE_FACTOR_HPP

#include <macros.hpp>
#include <table_function.hpp>
#include <cosmological_params.hpp>

class ScaleFactorFunctionClass : public Math::RealFunction
{
public:
	ScaleFactorFunctionClass() {}
	~ScaleFactorFunctionClass() {}
	
	void initialize(const CosmologicalParams& params);
	
	double age() const { check(!tableFunction_.empty(), ""); return tCut_ + tRadDom_; }
	
	double evaluate(double t) const;
	
	double time(double a) const; //Determine the time at which scale factor = a
	
	double comovingDistance(double t1, double t2) const;
	
	double hubble(double a) const;
    
    double growthFactor(double a) const;
	
private:
	typedef Math::TableFunction<double, double> TableFunctionType;
    
private:
    double growthFactorUnnormalizedRadDom(double a) const;
    
    double growthFactorFunction(double a) const;
	
private:
	double omegaM_;
	double omegaLambda_;
	double omegaR_;
	double omegaK_;
	double h0_;
	TableFunctionType tableFunction_;
    TableFunctionType growthFactor_;
    double growthFactorCoefficient_;
	double tCut_;
	double aCut_;
	double radDomCoeff_; // a = radDomCoeff_ * t^{1/2} in radiation domination epoch
	double tRadDom_;
};

#endif


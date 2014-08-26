#include <cmath>
#include <algorithm>
#include <set>

#include <unit_conversions.hpp>
#include <table_function.hpp>
#include <scale_factor.hpp>

using namespace Phys;
using namespace Math;

void
ScaleFactorFunctionClass::initialize(const CosmologicalParams& params)
{
	omegaM_ = params.getOmM();
	omegaLambda_ = params.getOmLambda();
	omegaR_ = params.getOmR();
	omegaK_ = params.getOmK();
	h0_ = params.getHubbleUnitless();

	const double deltaAFactor = double(1) / 10000;
	const double aFactor = 1 - deltaAFactor;
	const double aFactorSq = aFactor * aFactor;
	const double aFactorCube = aFactorSq * aFactor;
	const double aFactorForth = aFactorCube * aFactor;
	const double radiationDominationFactor = 10000;
	
	double a = 1;
	double t = 0;
	
	output_screen2("Calculating the scale factor..." << std::endl);
    
    std::set<double> aSet;
	
	TableFunctionType::iterator it = tableFunction_.insert(tableFunction_.end(), std::make_pair(t, a));
    aSet.insert(a);
    
	while(omegaR_ / omegaM_ < radiationDominationFactor)
	{
		double aPrime = a * h0_ * std::sqrt(omegaLambda_ + omegaM_ + omegaR_ + omegaK_);
		double deltaA = a * deltaAFactor;
		double deltaT = deltaA / aPrime;
		t += deltaT;
		tableFunction_.insert(it, std::make_pair(t, a));
        aSet.insert(a);
		
		a *= aFactor;
		omegaM_ /= aFactorCube;
		omegaR_ /= aFactorForth;
		omegaK_ /= aFactorSq;
	}
	
	output_screen2("Scale factor determined at " << tableFunction_.size() << " points!" << std::endl);
	
	tCut_ = t;
	aCut_ = a;
	
	radDomCoeff_ = std::sqrt(2 * h0_ * aCut_ * aCut_ * std::sqrt(omegaR_));
	tRadDom_ = 1 / (2 * h0_ * std::sqrt(omegaR_));

	omegaM_ = params.getOmM();
	omegaLambda_ = params.getOmLambda();
	omegaR_ = params.getOmR();
	omegaK_ = params.getOmK();
	h0_ = params.getHubbleUnitless();
    
    double gf = growthFactorUnnormalizedRadDom(aCut_);
    growthFactor_[aCut_] = gf;
    double aPrev = aCut_;
    double fPrev = growthFactorFunction(aPrev);
    for(std::set<double>::const_iterator it = aSet.begin(); it != aSet.end(); ++it)
    {
        double aNew = *it;
        double fNew = growthFactorFunction(aNew);
        gf += (fNew + fPrev) * (aNew - aPrev) / 2;
        growthFactor_[aNew] = gf;
        aPrev = aNew;
        fPrev = fNew;
    }
    
    growthFactorCoefficient_ = 1 / gf;
}

double
ScaleFactorFunctionClass::growthFactorFunction(double a) const
{
    double result = a * hubble(a) / h0_;
    result = 1 / (result * result * result);
    return result;
}

double
ScaleFactorFunctionClass::evaluate(double t) const
{
	check(!tableFunction_.empty(), "");
	
	if(t <= tRadDom_)
		return radDomCoeff_ * std::sqrt(t);
	
	t = tRadDom_ + tCut_ - t;
	check(t >= 0, "cannot evaluate at times bigger than the current age");
	check(t <= tCut_, "");
	
	return tableFunction_.evaluate(t);
}

double
ScaleFactorFunctionClass::time(double a) const
{
	check(!tableFunction_.empty(), "");
	check(a <= 1 && a >= 0, "invalid value for a");
	
	if(a <= aCut_)
		return a * a / (radDomCoeff_ * radDomCoeff_);
	
	double t1 = 0, t2 = tCut_;
	double tMiddle = (t1 + t2) / 2;
	double currentA = tableFunction_.evaluate(tMiddle);
	while(std::abs(a - currentA) / a > 0.0001)
	{
		if(a < currentA)
			t1 = tMiddle;
		else
			t2 = tMiddle;
		
		tMiddle = (t1 + t2) / 2;

		currentA = tableFunction_.evaluate(tMiddle);
	}
	return tRadDom_ + tCut_ - tMiddle;	
}

double
ScaleFactorFunctionClass::comovingDistance(double t1, double t2) const
{
	check(!tableFunction_.empty(), "");
	
	check(t1 < t2, "");
	check(t1 >= 0, "");
	check(t2 <= age(), t2 << ' ' << age());
	
	if(t2 <= tRadDom_)
		return 2 / radDomCoeff_ * (std::sqrt(t2) - std::sqrt(t1));
	
	double res = 0;
	if(t1 < tRadDom_)
	{
		res = 2 / radDomCoeff_ * (std::sqrt(tRadDom_) - std::sqrt(t1));
		t1 = tRadDom_;
	}
	
	t2 = age() - t2;
	t1 = age() - t1;
	
	std::swap(t1, t2);
	
	if(t1 < 0)
		t1 = 0;
	if(t2 > tCut_)
		t2 = tCut_;
	
	check(t1 <= t2, "");

	TableFunctionType::const_iterator beg = tableFunction_.lower_bound(t1), end = tableFunction_.lower_bound(t2);
	check(beg != tableFunction_.end(), "");
	check(end != tableFunction_.end(), "");
	if(beg == end)
		return res;
	
	double f1 = 1 / (*beg).second;
	double x1 = (*beg).first;
	
	while(beg++ != end)
	{
		double f2 = 1 / (*beg).second;
		double x2 = (*beg).first;
		res += (f1 + f2) / 2 * (x2 - x1);
		f1 = f2;
		x1 = x2;
	}
	
	return res;
}

double
ScaleFactorFunctionClass::hubble(double a) const
{
	check(!tableFunction_.empty(), "");
	
	check(a >= 0 && a <= 1, "invalid value");
	
	return h0_ * std::sqrt(omegaLambda_ + omegaM_ / (a * a * a) + omegaR_ / (a * a * a * a) + omegaK_ / (a * a));
}

double
ScaleFactorFunctionClass::growthFactorUnnormalizedRadDom(double a) const
{
    check(a >= 0 && a <= aCut_, "");
    double rdCoeff6 = radDomCoeff_ * radDomCoeff_ * radDomCoeff_;
    rdCoeff6 = rdCoeff6 * rdCoeff6;
    return 2 * a * a * a * a * h0_ * h0_ * h0_ / rdCoeff6;
}

double
ScaleFactorFunctionClass::growthFactor(double a) const
{
    check(a >= 0 && a <= 1, "");
    
    double result;
    if(a <= aCut_)
        result = growthFactorUnnormalizedRadDom(a);
    else
        result = growthFactor_.evaluate(a);
    
    result *= (growthFactorCoefficient_ * hubble(a) / h0_);
    return result;
}


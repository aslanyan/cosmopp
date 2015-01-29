#ifndef COSMO_PP_MARKOV_CHAIN_HPP
#define COSMO_PP_MARKOV_CHAIN_HPP

#include <vector>
#include <limits>
#include <ctime>

#include <function.hpp>
#include <table_function.hpp>
#include <random.hpp>

/// Posterior distribution for one parameter.
class Posterior1D : public Math::RealFunction
{
public:
    enum SmoothingMethod { GAUSSIAN_SMOOTHING = 0, SPLINE_SMOOTHING, SMOOTHING_MAX };
public:
    ///Constructior.
    Posterior1D(int seed = 0) : smooth_(NULL), cumulInv_(NULL), min_(std::numeric_limits<double>::max()), max_(-std::numeric_limits<double>::max()), minLike_(std::numeric_limits<double>::max()), generator_((seed == 0 ? std::time(0) : seed), 1e-5, 1.0 - 1e-5) {}

    /// Destructor.
    ~Posterior1D() { if(smooth_) delete smooth_; if(cumulInv_) delete cumulInv_; }

    /// Add a sample point. All of the sample points must be added before generating the distribution with generate.
    /// \param x The value of the parameter.
    /// \param prob The probability (weight) of the sample point.
    /// \param like The likelihood value for that given parameter value.
    void addPoint(double x, double prob, double like);

    /// Generate the distribution. This function should be called after all of the sample points have been added with addPoint.
    /// \param method The smoothing method. Can be GAUSSIAN_SMOOTHING for Gaussian smoothing or SPLINE_SMOOTHING for cubic spline smoothing.
    /// \param scale The smoothing scale. For Gaussian smoothing this is simply the smoothing scale. For spline smoothing this determines the distance between the points used for constructing the cubic spline. If not specified, the scale will be automatically determined from the number of sample points and their overall range.
    void generate(SmoothingMethod method = GAUSSIAN_SMOOTHING, double scale = 0);

    /// The minimum value of the parameter.
    double min() const { return min_; }

    /// The maximum value of the parameter.
    double max() const { return max_; }

    /// The value of the parameter for which the likelihood is maximum.
    double maxLikePoint() const { return maxLikePoint_; }
    
    /// The mean value of the parameter. Must be called after calling generate.
    double mean() const { check(smooth_, "not generated"); return mean_; }

    /// The value of the parameter where the distribution reaches its peak. Must be called after calling generate.
    double peak() const;

    /// The median value of the parameter. This is determined from the smoothed distribution. Must be called after calling generate.
    double median() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate(0.5 * norm_); }

    /// Get the one sigma range of the parameter assuming a two-sided distribution. This will be a one sigma range around the median. Must be called after calling generate.
    /// \param lower Will contain the lower point of the range upon return.
    /// \param upper Will contain the upper point of the range upon return.
    void get1SigmaTwoSided(double& lower, double& upper) const { check(cumulInv_, "not generated"); upper = cumulInv_->evaluate((1.0 - (1.0 - 0.683) / 2) * norm_); lower = cumulInv_->evaluate((1.0 - 0.683) / 2 * norm_); }

    /// Get the two sigma range of the parameter assuming a two-sided distribution. This will be a two sigma range around the median. Must be called after calling generate.
    /// \param lower Will contain the lower point of the range upon return.
    /// \param upper Will contain the upper point of the range upon return.
    void get2SigmaTwoSided(double& lower, double& upper) const { check(cumulInv_, "not generated"); upper = cumulInv_->evaluate((1.0 - (1.0 - 0.955) / 2) * norm_); lower = cumulInv_->evaluate((1.0 - 0.955) / 2 * norm_); }

    /// One sigma upper bound of the parameter. Must be called after calling generate.
    double get1SigmaUpper() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate(0.683 * norm_); }

    /// Two sigma upper bound of the parameter. Must be called after calling generate.
    double get2SigmaUpper() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate(0.955 * norm_); }

    /// One sigma lower bound of the parameter. Must be called after calling generate.
    double get1SigmaLower() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate((1.0 - 0.683) * norm_); }

    /// Two sigma lower bound of the parameter. Must be called after calling generate.
    double get2SigmaLower() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate((1.0 - 0.955) * norm_); }

    /// Calculate the distribution at a given point. Must be called after calling generate.
    virtual double evaluate(double x) const { check(smooth_, "not generated"); const double res = smooth_->evaluate(x) / norm_; return (res >= 0.0 ? res : 0.0); }

    /// Generate a random sample from this distribution
    double generateSample() { return cumulInv_->evaluate(generator_.generate() * norm_); }

    /// Write the distribution into a text file.
    /// \param fileName The name of the file.
    /// \param n The number of points (10,000 by default).
    void writeIntoFile(const char* fileName, int n = 10000) const;

private:
    double min_, max_;
    double minLike_, maxLikePoint_;
    double mean_;
    std::vector<double> points_, probs_;
    Math::RealFunction* smooth_;
    Math::TableFunction<double, double>* cumulInv_;
    double norm_;

    Math::UniformRealGenerator generator_;
};

/// Posterior distribution for two parameters. This is useful for making contour plots.
class Posterior2D : public Math::Function2<double, double, double>
{
public:
    /// Constructor.
    Posterior2D() : smooth_(NULL), cumulInv_(NULL), min1_(std::numeric_limits<double>::max()), min2_(std::numeric_limits<double>::max()), max1_(-std::numeric_limits<double>::max()), max2_(-std::numeric_limits<double>::max()), minLike_(std::numeric_limits<double>::max()) {}

    /// Destructor.
    ~Posterior2D() { if(smooth_) delete smooth_; if(cumulInv_) delete cumulInv_; }

    /// Add a sample point. All of the sample points must be added before generating the distribution with generate.
    /// \param x1 The value of the first parameter.
    /// \param x2 The value of the second parameter.
    /// \param prob The probability (weight) of the sample point.
    /// \param like The likelihood value for those given parameter values.
    void addPoint(double x1, double x2, double prob, double like);

    /// Generate the distribution. This function should be called after all of the sample points have been added with addPoint. The distribution is smoothed using two dimensional Gaussian smoothing.
    /// \param scale1 The smoothing scale for parameter 1. If not specified, the scale will be automatically determined from the number of sample points and their overall range.
    /// \param scale2 The smoothing scale for parameter 2. If not specified, the scale will be automatically determined from the number of sample points and their overall range.
    void generate(double scale1 = 0, double scale2 = 0);

    /// The minimum value of the first parameter.
    double min1() const { return min1_; }

    /// The minimum value of the second parameter.
    double min2() const { return min2_; }

    /// The maximum value of the first parameter.
    double max1() const { return max1_; }

    /// The maximum value of the second parameter.
    double max2() const { return max2_; }

    /// The value of the parameters for which the likelihood is maximum.
    void maxLikePoint(double& x1, double& x2) const { x1 = maxLikePoint1_; x2 = maxLikePoint2_; }

    /// Calculate the distribution at a given point. Must be called after calling generate.
    virtual double evaluate(double x1, double x2) const { check(smooth_, "not generated"); const double res = smooth_->evaluate(x1, x2) / norm_; return (res >= 0.0 ? res : 0.0); }

    /// The value of the distribution at the boundary of a region with a given confidence level.
    double getLevel(double confidence) const { check(cumulInv_, "not generated"); check (confidence >= 0 && confidence <= 1, "invalid confidence " << confidence); return cumulInv_->evaluate(confidence); }

    /// The value of the distribution at the boundary of the one sigma confidence region.
    double get1SigmaLevel() const { return getLevel(0.683); }

    /// The value of the distribution at the boundary of the two sigma confidence region.
    double get2SigmaLevel() const { return getLevel(0.955); }

    /// Write the distribution into a text file.
    /// \param fileName The name of the file.
    /// \param n The number of points in each dimension (1000 by default).
    void writeIntoFile(const char* fileName, int n = 1000) const;

private:
    double min1_, min2_, max1_, max2_;
    double minLike_, maxLikePoint1_, maxLikePoint2_;
    std::vector<double> points1_, points2_, probs_;
    Math::Function2<double, double, double>* smooth_;
    Math::TableFunction<double, double>* cumulInv_;
    double norm_;
};

/// A class for analyzing one or more Markov chains.
class MarkovChain
{
public:
    struct Element
    {
        double prob;
        double like;
        std::vector<double> params;

        bool operator <(const Element& other) const
        {
            return like < other.like;
        }
    };

public:
    /// Constructor for the case of a single chain.
    /// \param fileName The name of the file containing the chain.
    /// \param burnin The number of elements to ignore from the beginning of the chain.
    /// \param thin The thinning factor. Must be positive.
    MarkovChain(const char* fileName, unsigned long burnin = 0, unsigned int thin = 1);

    /// Constructor for the case of multiple chains.
    /// \param nChains The number of chains.
    /// \param fileNameRoot The root of the names of the files containing the chains. The actual file names should be this root followed by _ then the index of the chain (from 0 to nChains - 1) and then .txt
    /// \param burnin The number of elements to ignore from the beginning of the chain.
    /// \param thin The thinning factor. Must be positive.
    MarkovChain(int nChains, const char* fileNameRoot, unsigned long burnin = 0, unsigned int thin = 1);

    /// Destructor.
    ~MarkovChain();

    /// Allows to add another chain file.
    /// \param fileName The name of the file containing the chain.
    /// \param burnin The number of elements to ignore from the beginning of the chain.
    /// \param thin The thinning factor. Must be positive.
    void addFile(const char* fileName, unsigned long burnin = 0, unsigned int thin = 1);

    /// Returns the number of parameters.
    int nParams() const { return nParams_; }
    
    /// Get the one dimensional marginalized posterior distribution for a given parameter.
    /// \param paramIndex The index of the parameter, starting from 0.
    /// \param method The smoothing method. Can be Posterior1D::GAUSSIAN_SMOOTHING for Gaussian smoothing or Posterior1D::SPLINE_SMOOTHING for cubic spline smoothing.
    /// \param scale The smoothing scale. For Gaussian smoothing this is simply the smoothing scale. For spline smoothing this determines the distance between the points used for constructing the cubic spline. If not specified, the scale will be automatically determined from the number of sample points and their overall range.
    /// \return A pointer to the generated distribution. Must be deleted after using.
    Posterior1D* posterior(int paramIndex, Posterior1D::SmoothingMethod method = Posterior1D::GAUSSIAN_SMOOTHING, double scale = 0) const;

    /// Get the two dimensional marginalized posterior distribution for given parameters. This can be used to make contour plots. The distribution will be smoothed with Gaussian smoothing.
    /// \param paramIndex1 The index of the first parameter, starting from 0.
    /// \param paramIndex2 The index of the second parameter, starting from 0.
    /// \param scale1 The smoothing scale for parameter 1. If not specified, the scale will be automatically determined from the number of sample points and their overall range.
    /// \param scale2 The smoothing scale for parameter 2. If not specified, the scale will be automatically determined from the number of sample points and their overall range.
    /// \return A pointer to the generated distribution. Must be deleted after using.
    Posterior2D* posterior(int paramIndex1, int paramIndex2, double scale1 = 0, double scale2 = 0) const;

    /// -2ln(likelihood) for the maximum likelihood point.
    double maxLike() const { return minLike_; }

    /// Get a range of the points from at a given confidence level. For example, to get the one sigma range use pUpper = 0.683, pLower = 0 (the default values).
    /// \param container A vector where the elements will be written.
    /// \param pUpper The upper end of the confidence range.
    /// \param pLower The lower end of the confidence range.
    void getRange(std::vector<Element*>& container, double pUpper = 0.683, double pLower = 0) const;

private:
    void readFile(const char* fileName, unsigned long burnin, unsigned int thin, std::vector<Element*>& bigChain, double& maxP);
    void filterChain(std::vector<Element*>& bigChain, double minP);

private:
    std::vector<Element*> chain_;
    int nParams_;
    double minLike_;
};

#endif


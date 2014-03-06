#ifndef COSMO_PP_MARKOV_CHAIN_HPP
#define COSMO_PP_MARKOV_CHAIN_HPP

#include <vector>
#include <limits>

#include <function.hpp>
#include <table_function.hpp>

class Posterior1D : public Math::RealFunction
{
public:
    enum SmoothingMethod { GAUSSIAN_SMOOTHING = 0, SPLINE_SMOOTHING, SMOOTHING_MAX };
public:
    Posterior1D() : smooth_(NULL), cumul_(NULL), cumulInv_(NULL), min_(std::numeric_limits<double>::max()), max_(std::numeric_limits<double>::min()), minLike_(std::numeric_limits<double>::max()) {}
    ~Posterior1D() { if(smooth_) delete smooth_; if(cumul_) delete cumul_; if(cumulInv_) delete cumulInv_; }

    void addPoint(double x, double prob, double like);
    void generate(SmoothingMethod method = GAUSSIAN_SMOOTHING, double scale = 0);

    double min() const { return min_; }
    double max() const { return max_; }
    double maxLikePoint() const { return maxLikePoint_; }
    double mean() const { check(smooth_, "not generated"); return mean_; }
    double median() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate(0.5 * norm_); }
    void get1SigmaTwoSided(double& lower, double& upper) const { check(cumulInv_, "not generated"); upper = cumulInv_->evaluate((1.0 - (1.0 - 0.683) / 2) * norm_); lower = cumulInv_->evaluate((1.0 - 0.683) / 2 * norm_); }
    void get2SigmaTwoSided(double& lower, double& upper) const { check(cumulInv_, "not generated"); upper = cumulInv_->evaluate((1.0 - (1.0 - 0.955) / 2) * norm_); lower = cumulInv_->evaluate((1.0 - 0.955) / 2 * norm_); }
    double get1SigmaUpper() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate(0.683 * norm_); }
    double get2SigmaUpper() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate(0.955 * norm_); }
    double get1SigmaLower() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate((1.0 - 0.683) * norm_); }
    double get2SigmaLower() const { check(cumulInv_, "not generated"); return cumulInv_->evaluate((1.0 - 0.955) * norm_); }

    virtual double evaluate(double x) const { check(smooth_, "not generated"); const double res = smooth_->evaluate(x) / norm_; return (res >= 0.0 ? res : 0.0); }

private:
    double min_, max_;
    double minLike_, maxLikePoint_;
    double mean_;
    std::vector<double> points_, probs_, likes_;
    Math::RealFunction* smooth_;
    Math::TableFunction<double, double>* cumul_, *cumulInv_;
    double norm_;
};

class Posterior2D : public Math::Function2<double, double, double>
{
public:
    Posterior2D() : smooth_(NULL), cumul_(NULL), cumulInv_(NULL), min1_(std::numeric_limits<double>::max()), min2_(std::numeric_limits<double>::max()), max1_(std::numeric_limits<double>::min()), max2_(std::numeric_limits<double>::min()), minLike_(std::numeric_limits<double>::max()) {}
    ~Posterior2D() { if(smooth_) delete smooth_; if(cumul_) delete cumul_; if(cumulInv_) delete cumulInv_; }

    void addPoint(double x1, double x2, double prob, double like);
    void generate(double scale1 = 0, double scale2 = 0);

    double min1() const { return min1_; }
    double min2() const { return min2_; }
    double max1() const { return max1_; }
    double max2() const { return max2_; }
    void maxLikePoint(double& x1, double& x2) const { x1 = maxLikePoint1_; x2 = maxLikePoint2_; }

    virtual double evaluate(double x1, double x2) const { check(smooth_, "not generated"); const double res = smooth_->evaluate(x1, x2) / norm_; return (res >= 0.0 ? res : 0.0); }

    double getLevel(double confidence) const { check(cumulInv_, "not generated"); check (confidence >= 0 && confidence <= 1, "invalid confidence " << confidence); return cumulInv_->evaluate(confidence); }
    double get1SigmaLevel() const { return getLevel(0.683); }
    double get2SigmaLevel() const { return getLevel(0.955); }

private:
    double min1_, min2_, max1_, max2_;
    double minLike_, maxLikePoint1_, maxLikePoint2_;
    std::vector<double> points1_, points2_, probs_, likes_;
    Math::Function2<double, double, double>* smooth_;
    Math::TableFunction<double, double>* cumul_, *cumulInv_;
    double norm_;
};

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
    MarkovChain(const char* fileName);
    ~MarkovChain();

    int nParams() const { return nParams_; }
    Posterior1D* posterior(int paramIndex, Posterior1D::SmoothingMethod method = Posterior1D::GAUSSIAN_SMOOTHING, double scale = 0) const;
    Posterior2D* posterior(int paramIndex1, int paramIndex2, double scale1 = 0, double scale2 = 0) const;

    void getRange(std::vector<Element*>& container, double pUpper = 0.683, double pLower = 0) const;
private:
    std::vector<Element*> chain_;
    int nParams_;
};

#endif


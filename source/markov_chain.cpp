#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cubic_spline.hpp>
#include <gauss_smooth.hpp>
#include <progress_meter.hpp>
#include <markov_chain.hpp>

void
Posterior1D::addPoint(double x, double prob, double like)
{
    points_.push_back(x);
    probs_.push_back(prob);
    likes_.push_back(like);

    if(x < min_)
        min_ = x;
    if(x > max_)
        max_ = x;

    if(like < minLike_)
    {
        minLike_ = like;
        maxLikePoint_ = x;
    }
}

void
Posterior1D::generate(SmoothingMethod method, double scale)
{
    check(points_.size() >= 2, "at least 2 different data points need to be added before generating");
    check(points_.size() == probs_.size() && probs_.size() == likes_.size(), "");
    check(max_ > min_, "at least 2 different data points need to be added before generating");

    check(scale >= 0, "invalid scale " << scale);
    int resolution = 5 * (int)std::floor(std::sqrt(double(points_.size())));
    if(scale != 0)
        resolution = 5 * (int)std::floor((max_ - min_) / scale);
    else
        scale = 5 * (max_ - min_) / resolution;

    check(resolution > 0, "");

    std::vector<double> x(resolution + 2), y(resolution + 2, 0);
    const double d = (max_ - min_) / resolution;
    x[0] = min_;
    x[resolution + 1] = max_;
    for(int i = 0; i < resolution; ++i)
        x[i + 1] = min_ + d * i + d / 2;

    mean_ = 0;
    double totalP = 0;
    for(unsigned long i = 0; i < points_.size(); ++i)
    {
        const double p = points_[i]; 
        check(p >= min_, "");
        int k = (int)std::floor((p - min_) / d);
        check(k >= 0, "");
        if(k >= resolution)
            k = resolution - 1;

        y[k + 1] += probs_[i];
        mean_ += p * probs_[i];
        totalP += probs_[i];
    }

    check(totalP, "");
    mean_ /= totalP;

    if(smooth_)
    {
        check(cumul_, "");
        check(cumulInv_, "");
        delete smooth_;
        delete cumul_;
        delete cumulInv_;
    }

    switch(method)
    {
    case GAUSSIAN_SMOOTHING:
        smooth_ = new Math::GaussSmooth(x, y, scale);
        break;
    case SPLINE_SMOOTHING:
        smooth_ = new Math::CubicSpline(x, y);
        break;
    default:
        check(false, "");
        break;
    }

    const int N = 100 * resolution;
    const double delta = (x[x.size() - 1] - x[0]) / N;
    cumul_ = new Math::TableFunction<double, double>;
    cumulInv_ = new Math::TableFunction<double, double>;
    double maxVal = 0, maxX;
    norm_ = 0;
    for(int i = 0; i <= N; ++i)
    {
        double v = x[0] + i * delta;
        
        if(i == N)
            v = x[x.size() - 1];

        check(v <= x[x.size() - 1], "");

        double y = smooth_->evaluate(v);
        if(y < 0)
            y = 0;

        if(y > maxVal)
        {
            maxVal = y;
            maxX = v;
        }

        norm_ += y * delta;
        (*cumul_)[v] = norm_;
        (*cumulInv_)[norm_] = v;
    }
}

double
Posterior1D::peak() const
{
    const int nPoints = 10000;
    const double delta = (max() - min()) / nPoints;
    double peakVal = 0, x = min();
    for(int i = 0; i <= nPoints; ++i)
    {
        double t = min() + i * delta;
        if(i == nPoints)
            t = max();
        const double y = evaluate(t);
        if(y > peakVal)
        {
            peakVal = y;
            x = t;
        }
    }

    return x;
}

void
Posterior2D::addPoint(double x1, double x2, double prob, double like)
{
    points1_.push_back(x1);
    points2_.push_back(x2);
    probs_.push_back(prob);
    likes_.push_back(like);

    if(x1 < min1_)
        min1_ = x1;
    if(x2 < min2_)
        min2_ = x2;
    if(x1 > max1_)
        max1_ = x1;
    if(x2 > max2_)
        max2_ = x2;

    if(like < minLike_)
    {
        minLike_ = like;
        maxLikePoint1_ = x1;
        maxLikePoint2_ = x2;
    }
}

void
Posterior2D::generate(double scale1, double scale2)
{
    const int res = 10 * (int)std::ceil(std::sqrt(std::sqrt(double(points1_.size()))));
    int res1 = res, res2 = res;
    if(scale1 == 0)
        scale1 = 5 * (max1_ - min1_) / res1;
    else
        res1 = 5 * (int)std::floor((max1_ - min1_) / scale1);

    if(scale2 == 0)
        scale2 = 5 * (max2_ - min2_) / res2;
    else
        res2 = 5 * (int)std::floor((max2_ - min2_) / scale2);

    check(res1 > 0, "");
    check(res2 > 0, "");

    std::vector<double> x1(res1 + 2), x2(res2 + 2);
    std::vector<std::vector<double> > y(res1 + 2);
    for(int i = 0; i < res1 + 2; ++i)
        y[i].resize(res2 + 2, 0);

    const double d1 = (max1_ - min1_) / res1;
    x1[0] = min1_;
    x1[res1 + 1] = max1_;
    for(int i = 0; i < res1; ++i)
        x1[i + 1] = min1_ + d1 * i + d1 / 2;

    const double d2 = (max2_ - min2_) / res2;
    x2[0] = min2_;
    x2[res2 + 1] = max2_;
    for(int i = 0; i < res2; ++i)
        x2[i + 1] = min2_ + d2 * i + d2 / 2;

    for(unsigned long i = 0; i < points1_.size(); ++i)
    {
        const double p1 = points1_[i];
        const double p2 = points2_[i];
        check(p1 >= min1_, "");
        check(p2 >= min2_, "");
        int k1 = (int)std::floor((p1 - min1_) / d1);
        int k2 = (int)std::floor((p2 - min2_) / d2);
        check(k1 >= 0, "");
        check(k2 >= 0, "");
        if(k1 >= res1)
            k1 = res1 - 1;
        if(k2 >= res2)
            k2 = res2 - 1;

        y[k1 + 1][k2 + 1] += probs_[i];
    }

    if(smooth_)
    {
        check(cumul_, "");
        check(cumulInv_, "");

        delete smooth_;
        delete cumul_;
        delete cumulInv_;
    }

    std::vector<double> pointsX1, pointsX2, pointsY;
    for(int i = 0; i < x1.size(); ++i)
    {
        for(int j = 0; j < x2.size(); ++j)
        {
            pointsX1.push_back(x1[i]);
            pointsX2.push_back(x2[j]);
            pointsY.push_back(y[i][j]);
        }
    }

    smooth_ = new Math::GaussSmooth2D(pointsX1, pointsX2, pointsY, scale1, scale2);
    const int N1 = 10 * res1, N2 = 10 * res2;
    const double delta1 = (x1[x1.size() - 1] - x1[0]) / N1;
    const double delta2 = (x2[x2.size() - 1] - x2[0]) / N2;

    std::vector<double> probs;
    norm_ = 0;
    output_screen("Sampling the 2D distribution..." << std::endl);
    ProgressMeter met((N1 + 1) * (N2 + 1));
    for(int i = 0; i <= N1; ++i)
    {
        double v1 = x1[0] + i * delta1;
        if(i == N1)
            v1 = x1[x1.size() - 1];

        check(v1 <= x1[x1.size() - 1], "");
        
        for(int j = 0; j <= N2; ++j)
        {
            double v2 = x2[0] + j * delta2;
            if(j == N2)
                v2 = x2[x2.size() - 1];

            check(v2 <= x2[x2.size() - 1], "");

            double y = smooth_->evaluate(v1, v2);
            probs.push_back(y);
            norm_ += y * delta1 * delta2;
            met.advance();
        }
    }
    output_screen("OK" << std::endl);

    check(norm_ > 0, "");

    std::sort(probs.begin(), probs.end());
    cumul_ = new Math::TableFunction<double, double>;
    cumulInv_ = new Math::TableFunction<double, double>;

    double total = 0;
    int index = probs.size() - 1;
    while(total < 1 && index >= 0)
    {
        const double currentP = probs[index] / norm_;
        (*cumul_)[currentP] = total;
        (*cumulInv_)[total] = currentP;
        total += currentP * delta1 * delta2;
        --index;
    }
    (*cumul_)[0] = 1;
    (*cumulInv_)[1] = 0;
}

bool lessMarkovChainElementPointer(MarkovChain::Element* i, MarkovChain::Element* j)
{
    return (*i) < (*j);
}

MarkovChain::MarkovChain(const char* fileName, unsigned long burnin)
{
    minLike_ = std::numeric_limits<double>::max();
    addFile(fileName, burnin);
}

MarkovChain::MarkovChain(int nChains, const char* fileNameRoot, unsigned long burnin)
{
    check(nChains > 0, "need at least 1 chain");

    minLike_ = std::numeric_limits<double>::max();

    std::vector<Element*> bigChain;
    double maxP = std::numeric_limits<double>::min();
    for(int i = 0; i < nChains; ++i)
    {
        std::stringstream fileName;
        fileName << fileNameRoot;
        if(nChains > 1)
            fileName << '_' << i;
        fileName << ".txt";
        double thisMaxP;
        readFile(fileName.str().c_str(), burnin, bigChain, thisMaxP);
        if(thisMaxP > maxP)
            maxP = thisMaxP;
    }

    double minP = maxP / bigChain.size() / 1000;
    filterChain(bigChain, minP);

    output_screen("Sorting the chain..." << std::endl);
    std::sort(chain_.begin(), chain_.end(), lessMarkovChainElementPointer);
    output_screen("OK" << std::endl);
}

MarkovChain::~MarkovChain()
{
    for(unsigned long i = 0; i < chain_.size(); ++i)
        delete chain_[i];
}

void
MarkovChain::addFile(const char* fileName, unsigned long burnin)
{
    std::vector<Element*> bigChain;
    double maxP;
    readFile(fileName, burnin, bigChain, maxP);
    double minP = maxP / bigChain.size() / 1000;
    filterChain(bigChain, minP);

    output_screen("Sorting the chain..." << std::endl);
    std::sort(chain_.begin(), chain_.end(), lessMarkovChainElementPointer);
    output_screen("OK" << std::endl);
}

void
MarkovChain::readFile(const char* fileName, unsigned long burnin, std::vector<Element*>& bigChain, double& maxP)
{
    StandardException exc;
    std::ifstream in(fileName);

    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    output_screen("Reading the chain from file " << fileName << "..." << std::endl);
    nParams_ = -1;
    unsigned long line = 0;
    maxP = std::numeric_limits<double>::min();

    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);
        if(s == "")
            break;

        std::stringstream str(s);
        Element* elem = new Element;
        str >> elem->prob >> elem->like;

        if(elem->prob > maxP)
            maxP = elem->prob;

        if(elem->like < minLike_)
            minLike_ = elem->like;

        while(!str.eof())
        {
            double val = std::numeric_limits<double>::min();
            str >> val;
            if(val == std::numeric_limits<double>::min())
                break;

            elem->params.push_back(val);
        }

        if(nParams_ == -1)
            nParams_ = elem->params.size();

        else if (nParams_ != elem->params.size())
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid chain file " << fileName << ". There are " << elem->params.size() << " parameters on line " << line << " while the previous lines had " << nParams_ << " parameters.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(line >= burnin)
            bigChain.push_back(elem);

        ++line;
    }
    output_screen("OK" << std::endl);
    output_screen("Successfully read the chain. It has " << bigChain.size() << " elements, " << nParams_ << " parameters." << std::endl);
}

void
MarkovChain::filterChain(std::vector<Element*>& bigChain, double minP)
{
    output_screen("Filtering the chain..." << std::endl);
    ProgressMeter meter(bigChain.size());
    unsigned long left = 0;
    for(unsigned long i = 0; i < bigChain.size(); ++i)
    {
        Element* e = bigChain[i];
        if(e->prob < minP)
            delete e;
        else
        {
            chain_.push_back(e);
            ++left;
        }

        meter.advance();
    }
    output_screen("OK" << std::endl);
    output_screen(left << " elements left after filtering!" << std::endl);
}

Posterior1D*
MarkovChain::posterior(int paramIndex, Posterior1D::SmoothingMethod method, double scale) const
{
    check(paramIndex >= 0 && paramIndex < nParams_, "invalid parameter index " << paramIndex);
    Posterior1D* post = new Posterior1D;

    for(unsigned long i = 0; i < chain_.size(); ++i)
        post->addPoint(chain_[i]->params[paramIndex], chain_[i]->prob, chain_[i]->like);

    post->generate(method, scale);
    return post;
}

Posterior2D*
MarkovChain::posterior(int paramIndex1, int paramIndex2, double scale1, double scale2) const
{
    check(paramIndex1 >= 0 && paramIndex1 < nParams_, "invalid parameter index " << paramIndex1);
    check(paramIndex2 >= 0 && paramIndex2 < nParams_, "invalid parameter index " << paramIndex2);

    Posterior2D* post = new Posterior2D;
    for(unsigned long i = 0; i < chain_.size(); ++i)
        post->addPoint(chain_[i]->params[paramIndex1], chain_[i]->params[paramIndex2], chain_[i]->prob, chain_[i]->like);

    post->generate(scale1, scale2);
    return post;
}

void
MarkovChain::getRange(std::vector<Element*>& container, double pUpper, double pLower) const
{
    check(pUpper >= 0 && pUpper <= 1, "invalid probability " << pUpper << ", should be between 0 and 1");
    check(pLower >= 0 && pLower <= pUpper, "invalid lower probability " << pLower << ", should be between 0 and " << pUpper);
    container.clear();

    if(pUpper == 0)
        return;

    if(pLower == 1)
    {
        container.insert(container.end(), chain_.begin(), chain_.end());
        return;
    }

    double total = 0;
    std::vector<Element*>::const_iterator it = chain_.begin();
    while(total <= pUpper && it != chain_.end())
    {
        total += (*it)->prob;
        ++it;
        if(total > pLower)
            container.push_back(*it);
    }
}

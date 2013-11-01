#ifndef COSMO_CPP_RANDOM_HPP
#define COSMO_CPP_RANDOM_HPP

#include <numeric>

#include <boost/math/distributions/poisson.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>

namespace Math
{

/// Poisson distribution generator.
class PoissonGenerator
{
public:
    /// Constructor.
    /// \param seed The seed to use for the generator.
    /// \param mean The mean of the Poisson distribution.
    PoissonGenerator(int seed, double mean);

    /// Destructor.
    ~PoissonGenerator() {delete generator_;}
    
    /// A function to generate a random number from the distribution.
    /// \return A random number from the Poisson distribution.
    int generate();
    
private:
    std::vector<double> cumulative_;
    boost::mt19937 gen_;
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> >* generator_;
};

inline
PoissonGenerator::PoissonGenerator(int seed, double mean) : gen_(seed)
{
    using namespace boost::math;
    using namespace boost::random;
    
    std::vector<double> probabilities;
    poisson pd(mean);
    for(int i = 0; true; ++i)
    {
        double prob = pdf(pd, i);
        if(prob < 1e-10 && !probabilities.empty())
            break;
        
        probabilities.push_back(prob);
    }
    
    std::partial_sum(probabilities.begin(), probabilities.end(), std::back_inserter(cumulative_));
    boost::uniform_real<> dist(0, cumulative_.back());
    
    generator_ = new boost::variate_generator<boost::mt19937&, boost::uniform_real<> >(gen_, dist);
}

inline int
PoissonGenerator::generate()
{
    return (std::lower_bound(cumulative_.begin(), cumulative_.end(), (*generator_)()) - cumulative_.begin());
}
    
} // namespace Math

#endif

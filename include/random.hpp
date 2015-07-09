#ifndef COSMO_PP_RANDOM_HPP
#define COSMO_PP_RANDOM_HPP

#include <random>

namespace Math
{

/// Uniform int generator
class UniformIntGenerator
{
public:
    /// Constructor.
    /// \param seed The seed to use for the generator.
    /// \param min The minimum of the range.
    /// \param max The maximum of the range.
    UniformIntGenerator(int seed, int min = 0, int max = 1) : gen_(seed), dist_(min, max) {}

    /// Destructor.
    ~UniformIntGenerator() {}

    /// A function to generate a random number form the distribution.
    /// \return A random number from the uniform distribution.
    int generate() { return dist_(gen_); }

private:
    std::mt19937 gen_;
    std::uniform_int_distribution<int> dist_;
};

/// Uniform real generator
class UniformRealGenerator
{
public:
    /// Constructor.
    /// \param seed The seed to use for the generator.
    /// \param min The minimum of the range.
    /// \param max The maximum of the range.
    UniformRealGenerator(int seed, double min = 0.0, double max = 1.0) : gen_(seed), dist_(min, max) {}

    /// Destructor.
    ~UniformRealGenerator() {}

    /// A function to generate a random number form the distribution.
    /// \return A random number from the uniform distribution.
    double generate() { return dist_(gen_); }

private:
    std::mt19937 gen_;
    std::uniform_real_distribution<double> dist_;
};

/// Gaussian distribution generator
class GaussianGenerator
{
public:
    /// Constructor.
    /// \param seed The seed to use for the generator.
    /// \param mean The mean of the Gaussian.
    /// \param sigma The sigma of the Gaussian.
    GaussianGenerator(int seed, double mean, double sigma) : gen_(seed), dist_(mean, sigma) {}

    /// Destructor.
    ~GaussianGenerator() {}

    /// A function to generate a random number from the distribution.
    /// \return A random number from the Gaussian distribution.
    double generate() { return dist_(gen_); }

private:
    std::mt19937 gen_;
    std::normal_distribution<> dist_;
};

/// Poisson distribution generator.
class PoissonGenerator
{
public:
    /// Constructor.
    /// \param seed The seed to use for the generator.
    /// \param mean The mean of the Poisson distribution.
    PoissonGenerator(int seed, double mean) : gen_(seed), dist_(mean) {}

    /// Destructor.
    ~PoissonGenerator() {}
    
    /// A function to generate a random number from the distribution.
    /// \return A random number from the Poisson distribution.
    int generate() { return dist_(gen_); }
    
private:
    std::mt19937 gen_;
    std::poisson_distribution<int> dist_;
};

} // namespace Math

#endif

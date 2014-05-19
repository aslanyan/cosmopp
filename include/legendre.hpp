#ifndef COSMO_PP_LEGENDRE_HPP
#define COSMO_PP_LEGENDRE_HPP

#include <vector>

#include <polynomial.hpp>

namespace Math
{

class Legendre
{
public:
    Legendre()
    {
        Polynomial p0(1), p1(2);
        p0.parameter(0) = 1;
        p1.parameter(1) = 1;

        p_.push_back(p0);
        p_.push_back(p1);
    }

    ~Legendre() {}

    const Polynomial& get(unsigned int l)
    {
        while(l >= p_.size())
            add();

        return p_[l];
    }

private:
    void add()
    {
        int n = p_.size();
        check(n >= 2, "");

        Polynomial pNew = ((2 * n - 1) * p_[1] * p_[n - 1] - (n - 1) * p_[n - 2]) / n;
        p_.push_back(pNew);
    }
private:
    std::vector<Polynomial> p_;
};

} // namespace Math

#endif


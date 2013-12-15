#ifndef COSMO_PP_WIGNER_3J_HPP
#define COSMO_PP_WIGNER_3J_HPP

#include <cmath>
#include <set>
#include <algorithm>

#include <macros.hpp>

namespace Math
{
    
/* struct Treats
 {
    //void totalNum(long n);
    void process(int l1, int l2, int l3, double val);
 };
*/
    
/// Wigner 3j symbols calculator with 0 m-s.

/// This class recursively calculates the Wigner 3j symbols with 0 m-s and integer l-s. This should be used when one needs them all at once, 
/// that's where the speedup makes a huge difference. 
/// The template Treats needs to have a process(int l1, int l2, int l3, double val) function which is called for all of the non-zero values of 
/// the Wigner symbol up to a given maximum from the calculate function.
template<typename Treats>
class Wigner3JZeroM
{
public:
    /// Constructor.
    /// \param t The Treats object.
    Wigner3JZeroM(Treats& t) : t_(t) {}
    
    /// The main calculator.
    /// \param LMaxHalf The half of the maximum of L = l1 + l2 + l3. The process function for the Treats object will be called for all the non-zero symbols up to that maximum.
    void calculate(int LMaxHalf);
    
private:
    void process(int l1, int l2, int l3, double val);
    Treats& t_;
};
    
inline double wigner3jFactorForNext(int LHalf, int indexGoingUp, int indexGoingDown)
{
    return std::sqrt(double(2 * LHalf - 2 * indexGoingDown + 2) * (2 * LHalf - 2 * indexGoingDown + 1) / (double(2 * LHalf - 2 * indexGoingUp) * (2 * LHalf - 2 * indexGoingUp - 1))) * double(LHalf - indexGoingUp) / double(LHalf - indexGoingDown + 1);
}

template<typename Treats>
void
Wigner3JZeroM<Treats>::calculate(int LMaxHalf)
{
    /*long total = 0;
    for(int LHalf = 0; LHalf <= LMaxHalf; ++LHalf)
    {
        const int l1Max = 2 * LHalf / 3;
        for(int l1 = 0; l1 <= l1Max; ++l1)
        {
            const int l2Min = std::max(l1, LHalf - l1);
            const int l2Max = (2 * LHalf - l1) / 2;
            check(l2Max >= l2Min, "");
            for(int l2 = l2Min; l2 <= l2Max; ++l2)
            {
                const int l3 = 2 * LHalf - l1 - l2;
                check(l3 >= l2, "");
                check(l3 <= l1 + l2, "");
                check(l1 <= l2, "");
                ++total;
            }
        }
    }
    
    t_.totalNum(total);*/
    
    double prevFirst = 1;
    
    process(0, 0, 0, 1);
    
    for(int LHalf = 1; LHalf <= LMaxHalf; ++LHalf)
    {
        prevFirst *= -std::sqrt(double(2 * (LHalf - 1) + 1) / (2 * (LHalf - 1) + 3));
        
        double first = prevFirst;
        int prevL2Min = LHalf;
        
        const int l1Max = 2 * LHalf / 3;
        process(0, LHalf, LHalf, first);
        
        for(int l1 = 1; l1 <= l1Max; ++l1)
        {
            const int l2Min = std::max(l1, LHalf - l1);
            const int l2Max = (2 * LHalf - l1) / 2;
            check(l2Max >= l2Min, "");
            
            check(l2Min == prevL2Min || l2Min == prevL2Min - 1 || l2Min == prevL2Min + 1, "");
            
            int prevIndexThatDecreases = prevL2Min;
            
            if(l2Min == prevL2Min || l2Min == prevL2Min + 1)
            {
                prevIndexThatDecreases = 2 * LHalf - prevL2Min - (l1 - 1);
            }
            
            first *= wigner3jFactorForNext(LHalf, l1 - 1, prevIndexThatDecreases);
            
            if(l2Min == prevL2Min + 1)
                first *= wigner3jFactorForNext(LHalf, prevL2Min, prevIndexThatDecreases - 1);
            
            double current = first;
            
            process(l1, l2Min, 2 * LHalf - l1 - l2Min, current);
            
            prevL2Min = l2Min;
            
            for(int l2 = l2Min + 1; l2 <= l2Max; ++l2)
            {
                const int l3 = 2 * LHalf - l1 - l2;
                current *= wigner3jFactorForNext(LHalf, l2 - 1, l3 + 1);
                check(l3 >= l2, "");
                check(l3 <= l1 + l2, "");
                check(l1 <= l2, "");
                process(l1, l2, l3, current);
            }
        }
    }
}

struct IntegerTriplet
{
    int l1, l2, l3;
    bool operator < (const IntegerTriplet& other) const
    {
        if(l1 < other.l1)
            return true;
        if(l1 > other.l1)
            return false;
        if(l2 < other.l2)
            return true;
        if(l2 > other.l2)
            return false;
        return l3 < other.l3;
    }
};
    
template<typename Treats>
void
Wigner3JZeroM<Treats>::process(int l1, int l2, int l3, double val)
{
    std::set<IntegerTriplet> triplets;
    IntegerTriplet t;
    t.l1 = l1;
    t.l2 = l2;
    t.l3 = l3;
    triplets.insert(t);
    std::swap(t.l1, t.l2);
    triplets.insert(t);
    std::swap(t.l2, t.l3);
    triplets.insert(t);
    std::swap(t.l1, t.l3);
    triplets.insert(t);
    std::swap(t.l1, t.l2);
    triplets.insert(t);
    std::swap(t.l2, t.l3);
    triplets.insert(t);

    for(std::set<IntegerTriplet>::const_iterator it = triplets.begin(); it != triplets.end(); ++it)
        t_.process((*it).l1, (*it).l2, (*it).l3, val);
}

} //namespace Math

#endif

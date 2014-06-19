#include <fstream>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <numerics.hpp>
#include <three_rotation.hpp>
#include <whole_matrix.hpp>

/*#include "healpix_base.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "healpix_map.h"
#include "rotmatrix.h"
#include "alm_powspec_tools.h"
#include "chealpix.h"*/

WholeMatrix::WholeMatrix(int lMin, int lMax) : lMin_(lMin), lMax_(lMax)
{
    initialize();
}

void
WholeMatrix::initialize()
{
    check(lMin_ >= 0, "");
    check(lMax_ >= lMin_, "");
    
    data_.resize(lMax_ + 1);
    for(int l1 = lMin_; l1 <= lMax_; ++l1)
    {
        data_[l1].resize(2 * l1 + 1);
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            data_[l1][l1 + m1].resize(lMax_ + 1);
            for(int l = lMin_; l <= lMax_; ++l)
            {
                data_[l1][l1 + m1][l].resize(2 * l + 1, 0);
            }
        }
    }
}

bool
WholeMatrix::checkIndices(int l, int m) const
{
    if(l < lMin_ || l > lMax_)
        return false;
    
    if(m < -l || m > l)
        return false;
    
    return true;
}

double
WholeMatrix::element(int l1, int m1, int l, int m) const
{
    check(checkIndices(l1, m1), "invalid l1, m1, l1 = " << l1 << " m1 = " << m1);
    check(checkIndices(l, m), "invalid l, m, l = " << l << " m = " << m);
    return data_[l1][l1 + m1][l][l + m];
}

double&
WholeMatrix::element(int l1, int m1, int l, int m)
{
    check(checkIndices(l1, m1), "invalid l1, m1, l1 = " << l1 << " m1 = " << m1);
    check(checkIndices(l, m), "invalid l, m, l = " << l << " m = " << m);
    return data_[l1][l1 + m1][l][l + m];
}

void
WholeMatrix::readFromTextFile(const char* fileName)
{
    StandardException exc;
    
    std::ifstream in(fileName);
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    std::string str;
    std::getline(in, str);
    std::stringstream sstr(str);
    sstr >> lMin_ >> lMax_;
    
    data_.clear();
    initialize();
    
    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);
        if(s == "")
            break;
        
        std::stringstream str(s);
        int l, m, l1, m1;
        double val;
        str >> l >> m >> l1 >> m1 >> val;
        
        element(l1, m1, l, m) = val;
    }
    in.close();
}

void
WholeMatrix::writeIntoTextFile(const char* fileName) const
{
    StandardException exc;
    std::ofstream out(fileName);
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    out << lMin_ << ' ' << lMax_ << std::endl;
    for(int l = lMin_; l <= lMax_; ++l)
    {
        for(int l1 = lMin_; l1 <= lMax_; ++l1)
        {
            for(int m = -l; m <= l; ++m)
            {
                for(int m1 = -l1; m1 <= l1; ++m1)
                    out << l << ' ' << m << ' ' << l1 << ' ' << m1 << ' ' << element(l1, m1, l, m) << std::endl;
            }
        }
    }
    out.close();
}

void
WholeMatrix::readFromFile(const char* fileName)
{
    StandardException exc;
    std::ifstream in(fileName, std::ios::in | std::ios::binary);
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    in.read((char*)(&lMin_), sizeof(int));
    in.read((char*)(&lMax_), sizeof(int));
    
    data_.clear();
    
    initialize();
    
    for(int l1 = lMin_; l1 <= lMax_; ++l1)
    {
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            for(int l = lMin_; l <= lMax_; ++l)
            {
                check(data_[l1][l1 + m1][l].size() == 2 * l + 1, "");
                in.read((char*)(&(data_[l1][l1 + m1][l][0])), (2 * l + 1) * sizeof(double));
            }
        }
    }
    in.close();
}

void
WholeMatrix::writeIntoFile(const char* fileName) const
{
    std::ofstream out(fileName, std::ios::binary | std::ios::out);
    if(!out)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    out.write((char*)(&lMin_), sizeof(int));
    out.write((char*)(&lMax_), sizeof(int));
    
    for(int l1 = lMin_; l1 <= lMax_; ++l1)
    {
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            for(int l = lMin_; l <= lMax_; ++l)
            {
                check(data_[l1][l1 + m1][l].size() == 2 * l + 1, "");
                out.write((char*)(&(data_[l1][l1 + m1][l][0])), (2 * l + 1) * sizeof(double));
            }
        }
    }
    out.close();
}


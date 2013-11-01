#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <utils.hpp>
#include <c_matrix.hpp>

#include "chealpix.h"

CMatrix::CMatrix(int nPix) : nPix_(nPix)
{
    initialize();
}

void
CMatrix::initialize()
{
    check(nPix_ > 0, "the number of pixels must be positive.");
    
    matrix_.resize(nPix_ * (nPix_ + 1) / 2, 0);
}

int
CMatrix::getIndex(int i, int j) const
{
    check(i >= 0 && i < nPix_, "invalid index" << i);
    check(j >= 0 && j < nPix_, "invalid index" << j);
    
    if(i > j)
        std::swap(i, j);
    
    const int index = j * (j + 1) / 2 + i;
    check(index < matrix_.size(), "");
    return index;
}

void
CMatrix::readFromFile(const char* fileName)
{
    std::ifstream inC;
    inC.open(fileName, std::ios::in | std::ios::binary);
    
    if(!inC)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Covariance matrix file " << fileName << " cannot be read.";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    inC.read((char*)(&nPix_), sizeof(int));
    
    initialize();
    inC.read((char*)(&(matrix_[0])), matrix_.size() * sizeof(double));
    int commentSize;
    inC.read((char*)(&commentSize), sizeof(int));
    comment_.resize(commentSize);
    inC.read((char*)(&(comment_[0])), commentSize * sizeof(char));
    inC.close();
}

void
CMatrix::writeIntoFile(const char* fileName) const
{
    std::ofstream cOut(fileName, std::ios::binary | std::ios::out);
    if(!cOut)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    cOut.write((char*)(&nPix_), sizeof(int));
    cOut.write((char*)(&(matrix_[0])), matrix_.size() * sizeof(double));
    int commentSize = comment_.size();
    cOut.write((char*)(&commentSize), sizeof(int));
    cOut.write((char*)(&(comment_[0])), commentSize * sizeof(char));
    cOut.close();
}

void
CMatrix::writeIntoTextFile(const char* fileName) const
{
    std::ofstream out(fileName);
    if(!out)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    out << nPix_ << std::endl;
    out << comment_ << std::endl;
    
    for(int j = 0; j < nPix_; ++j)
    {
        for(int i = 0; i <= j; ++i)
        {
            out << i << '\t' << j << '\t' << element(i, j) << std::endl;
        }
    }
    out.close();
}

void
CMatrix::readFromTextFile(const char* fileName)
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
    in >> nPix_;
    matrix_.clear();
    initialize();
    
    std::getline(in, comment_);
    
    while(!in.eof())
    {
        std::string s;
        in >> s;
        if(s == "")
            break;
        
        std::stringstream str(s);
        int i, j;
        double val;
        str >> i >> j >> val;
        if(i < 0 || i >= nPix_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid index i = " << i << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }
        if(j < 0 || j >= nPix_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid index j = " << j << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }
        element(i, j) = val;
    }
}

void
CMatrix::maskMatrix(const char* maskFileName)
{
    StandardException exc;
    
    long nSide;
    std::vector<int> goodPixels;
    Utils::readMask(maskFileName, nSide, goodPixels);
    
    const int nPix = (int)nside2npix(nSide);
    
    if(nPix_ != nPix)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The covariance matrix has " << nPix_ << " pixels, while there are " << nPix << " pixels in the mask. They need to be the same.";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    maskMatrix(goodPixels);
}

void
CMatrix::maskMatrix(const std::vector<int>& goodPixels)
{
    const int goodPixelsSize = goodPixels.size();
    
    std::vector<double> newMatrix(goodPixelsSize * (goodPixelsSize + 1) / 2);
        
    for(int j = 0; j < goodPixelsSize; ++j)
    {
        for(int i = 0; i <= j; ++i)
        {
            const int j1 = goodPixels[j];
            const int i1 = goodPixels[i];
            const int index = getIndex(i, j);
            newMatrix[index] = element(i1, j1);
        }
    }
    nPix_ = goodPixelsSize;
    matrix_ = newMatrix;
}

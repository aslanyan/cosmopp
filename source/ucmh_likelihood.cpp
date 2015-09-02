#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <numerics.hpp>
#include <ucmh_likelihood.hpp>

UCMHLikelihood::UCMHLikelihood(const char* fileName, bool lateKineticDecoupling) : pkMax_(1.0), likeMax_(1e10)
{
    output_screen("Initializing the UCMH likelihoods from the file: " << fileName << "..." << std::endl);

    std::ifstream in(fileName);
    StandardException exc;
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    std::string s;
    std::getline(in, s);

    if(s[0] != 'k' && !(s[0] == ' ' && s[1] == 'k'))
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid input file format. The first line should start with \"k\".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    std::vector<double> cl;
    extractCLValues(s, cl);
    createClToLike();

    std::vector<double> likes(cl.size());
    for(int i = 0; i < cl.size(); ++i)
    {
        likes[i] = clToLike_.evaluate(cl[i]);
        output_screen1("cl = " << cl[i] << " --> like = " << likes[i] << std::endl);
    }

    int line = 1;

    while(!in.eof())
    {
        std::getline(in, s);
        if(s == "")
            break;
        ++line;

        std::stringstream str(s);
        double k;
        str >> k;
        
        if(k < 1e-3 || k > 1e10)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid input on line " << line << " of the input file " << fileName << ". The k value of " << k << " is invalid.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(lateKineticDecoupling && k > 750)
            continue;

        double pkMin = 10 * pkMax_;

        for(int i = 0; i < likes.size(); ++i)
        {
            double pk;
            str >> pk;
            if((pk < 1e-10 || pk > pkMax_) && pk != 0)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Invalid input on line " << line << " of the input file " << fileName << ". The pk value of " << pk << " is invalid.";
                exc.set(exceptionStr.str());
                throw exc;
            }

            if(pk != 0)
            {
                likes_[k][pk] = likes[i];
                if(pk < pkMin)
                    pkMin = pk;
            }
        }

        check(pkMin <= pkMax_, "");

        likes_[k][0.9 * pkMin] = 0;
        likes_[k][pkMax_] = likeMax_;
    }

    output_screen("OK" << std::endl);
}

double
UCMHLikelihood::calculate(const Math::RealFunction& ps) const
{
    double likeMax = 0;
    for(auto it = likes_.cbegin(); it != likes_.cend(); ++it)
    {
        const double k = it->first;
        const LikelihoodType& like = it->second;

        const double pk = ps.evaluate(k);
        double l;
        if(pk > pkMax_)
            l = likeMax_;
        else
        {
            const double pkMin = like.cbegin()->first;
            if(pk <= pkMin)
                l = 0;
            else
                l = like.evaluate(pk);
        }

        if(l > likeMax)
            likeMax = l;
    }

    return likeMax;
}

void
UCMHLikelihood::extractCLValues(const std::string& s, std::vector<double>& cl) const
{
    cl.clear();
    std::size_t pos = 0;

    StandardException exc;
    while(true)
    {
        pos = s.find('=', pos);
        if(pos == std::string::npos)
            break;

        ++pos;

        std::size_t space = s.find(' ', pos);
        if(space == std::string::npos)
            space = s.size();

        std::stringstream str;
        str << s.substr(pos, space - pos);

        double c;
        str >> c;

        if(c <= 0 || c >= 1)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid value of cl = " << c << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        cl.push_back(c);

        if(space == s.size())
            break;
        
        pos = space;
    }

    output_screen("The following cl values have been read:");
    for(int i = 0; i < cl.size(); ++i)
    {
        output_screen(" " << cl[i]);
    }

    output_screen(std::endl);
}

void
UCMHLikelihood::createClToLike()
{
    clToLike_.clear();
    const double norm = 1.0 / std::sqrt(2.0 * Math::pi);
    double cl = 0;
    double x = 0;
    double y = norm;

    clToLike_[0] = 0;

    const double xMax = 10;
    const int n = 100000;
    const double deltaX = xMax / n;

    for(int i = 1; i <= n; ++i)
    {
        const double xNew = i * deltaX;
        const double l = xNew * xNew;
        const double yNew = std::exp(-l / 2.0) * norm;

        cl += (xNew - x) * (yNew + y);

        if(cl >= 1)
        {
            check(Math::areEqual(cl, 1.0, 1e-4), "");
            cl = 1;
        }
        clToLike_[cl] = l;
        y = yNew;
        x = xNew;
    }

    check(Math::areEqual(cl, 1.0, 1e-4), "");
}


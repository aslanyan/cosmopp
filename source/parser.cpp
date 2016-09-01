#include <fstream>
#include <sstream>
#include <cstdlib>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <parser.hpp>
#include <cosmo_mpi.hpp>

void
Parser::readFile(const char *fileName)
{
    std::ifstream in(fileName);
    StandardException exc;
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);

        std::string sClean;
        for(auto it = s.begin(); it != s.end(); ++it)
        {
            const char c = *it;
            if(c == ' ' || c == '\t')
                continue;
            sClean.push_back(c);
        }

        if(sClean.empty())
            continue;
        
        if(sClean[0] == '#')
            continue;

        auto equalSign = sClean.find('=');
        if(equalSign == std::string::npos)
            continue;

        if(equalSign + 1 == sClean.size())
        {
            if(CosmoMPI::create().isMaster())
            {
                output_screen("Invalid string: " << s << std::endl << "\tThere is no value after the equal sign. IGNORING!" << std::endl);
            }
            continue;
        }

        const std::string key = sClean.substr(0, equalSign);
        const std::string val = sClean.substr(equalSign + 1, sClean.size() - equalSign - 1);

        (*this)[key] = val;
    }
}

int
Parser::getInt(const std::string& s) const
{
    auto it = find(s);
    StandardException exc;
    if(it == cend())
    {
        std::stringstream exceptionStr;
        exceptionStr << "The key " << s << " is not found!";
        exc.set(exceptionStr.str());
        throw exc;
    }
    return std::atoi((it->second).c_str());
}

int
Parser::getInt(const std::string& s, int def)
{
    if(find(s) == end())
    {
        std::stringstream str;
        str << def;
        (*this)[s] = str.str();
        return def;
    }
    return getInt(s);
}

double
Parser::getDouble(const std::string& s) const
{
    auto it = find(s);
    StandardException exc;
    if(it == cend())
    {
        std::stringstream exceptionStr;
        exceptionStr << "The key " << s << " is not found!";
        exc.set(exceptionStr.str());
        throw exc;
    }
    return std::atof((it->second).c_str());
}

double
Parser::getDouble(const std::string& s, double def)
{
    if(find(s) == end())
    {
        std::stringstream str;
        str << def;
        (*this)[s] = str.str();
        return def;
    }
    return getDouble(s);
}

std::string
Parser::getStr(const std::string& s) const
{
    auto it = find(s);
    StandardException exc;
    if(it == cend())
    {
        std::stringstream exceptionStr;
        exceptionStr << "The key " << s << " is not found!";
        exc.set(exceptionStr.str());
        throw exc;
    }
    return it->second;
}

std::string
Parser::getStr(const std::string& s, const std::string& def)
{
    if(find(s) == end())
    {
        (*this)[s] = def;
        return def;
    }
    return getStr(s);
}

bool
Parser::getBool(const std::string& s) const
{
    auto it = find(s);
    StandardException exc;
    if(it == cend())
    {
        std::stringstream exceptionStr;
        exceptionStr << "The key " << s << " is not found!";
        exc.set(exceptionStr.str());
        throw exc;
    }
    const std::string& r = it->second;
    if(r == "yes" || r == "YES" || r == "Yes" || r == "y" || r == "Y" || r == "true" || r == "TRUE" || r == "True" || r == "t" || r == "T" || r == ".true." || r == ".True." || r == ".TRUE." || r == "1")
        return true;
    if(r == "no" || r == "NO" || r == "No" || r == "n" || r == "N" || r == "false" || r == "FALSE" || r == "False" || r == "f" || r == "F" || r == ".false." || r == ".False." || r == ".FALSE." || r == "0")
        return false;
    
    std::stringstream exceptionStr;
    exceptionStr << "Cannot interpret boolean value of \"" << r << "\" of the key \"" << s << "\"!";
    exc.set(exceptionStr.str());
    throw exc;
}

bool
Parser::getBool(const std::string& s, bool def)
{
    if(find(s) == end())
    {
        (*this)[s] = (def ? "true" : "false");
        return def;
    }
    try {
        return getBool(s);
    } catch (std::exception& e)
    {
        return def;
    }
}

void
Parser::dump() const
{
    if(!CosmoMPI::create().isMaster())
        return;
    output_screen("All parameters in parser:" << std::endl);
    for(auto it = cbegin(); it != cend(); ++it)
    {
        output_screen('\t' << it->first << " = " << it->second << std::endl);
    }
    output_screen("DONE" << std::endl);
}


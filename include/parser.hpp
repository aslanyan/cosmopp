#ifndef COSMO_PP_PARSER_HPP
#define COSMO_PP_PARSER_HPP

#include <string>
#include <map>

class Parser : public std::map<std::string, std::string>
{
public:
    Parser() {}

    Parser(const char *fileName)
    {
        readFile(fileName);
    }

    ~Parser() {}

    void dump() const;

    void readFile(const char* fileName);

    int getInt(const std::string& s) const;
    int getInt(const std::string& s, int def);

    double getDouble(const std::string& s) const;
    double getDouble(const std::string& s, double def);

    std::string getStr(const std::string& s) const;
    std::string getStr(const std::string& s, const std::string& def);
};

#endif


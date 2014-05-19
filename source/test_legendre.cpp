#include <test_legendre.hpp>
#include <legendre.hpp>

std::string
TestLegendre::name() const
{
    return std::string("LEGENDRE POLYNOMIAL TESTER");
}

unsigned int
TestLegendre::numberOfSubtests() const
{
    return 3;
}

void
TestLegendre::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 3, "invalid index " << i);
    
    using namespace Math;
    Legendre legendre;
    Polynomial p(1);
    switch(i)
    {
    case 0:
        subTestName = std::string("p4");
        p = legendre.get(4);
        res = p.evaluate(2);
        expected = 55.375;
        break;
    case 1:
        subTestName = std::string("degree");
        p = legendre.get(1000);
        res = p.numberOfParams() - 1;
        expected = 1000;
        break;
    case 2:
        subTestName = std::string("p10");
        p = legendre.get(10);
        res = p.parameter(8);
        expected = -double(109395) / 256;
        break;
    default:
        check(false, "");
        break;
    }
}

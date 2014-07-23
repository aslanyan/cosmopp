#include <macros.hpp>
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
    return 16;
}

void
TestLegendre::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < numberOfSubtests(), "invalid index " << i);
    
    using namespace Math;
    Legendre legendre;
    AssociatedLegendre associatedLegendre;
    switch(i)
    {
    case 0:
        subTestName = std::string("p4");
        res = legendre.calculate(4, 2);
        expected = 55.375;
        break;
    case 1:
        subTestName = std::string("p10");
        res = legendre.calculate(10, 0.5);
        expected = -0.188228607177734375;
        break;
    case 2:
        subTestName = std::string("p64");
        res = legendre.calculate(64, 0.1);
        expected = 0.098026402863;
        break;
    case 3:
        subTestName = std::string("p1000");
        res = legendre.calculate(1000, -0.25);
        expected = 0.0023444296560;
        break;
    case 4:
        subTestName = std::string("p10000");
        res = legendre.calculate(10000, -0.5);
        expected = -0.006062503808317;
        break;
    case 5:
        subTestName = std::string("associated_0_0");
        res = associatedLegendre.calculate(0, 0, 0.5);
        expected = 1.0;
        break;
    case 6:
        subTestName = std::string("associated_5_5");
        res = associatedLegendre.calculate(5, 5, -0.5);
        expected = -460.3466287;
        break;
    case 7:
        subTestName = std::string("associated_10_7");
        res = associatedLegendre.calculate(10, 7, 0.1);
        expected = 1.55806753e6;
        break;
    case 8:
        subTestName = std::string("associated_10_9");
        res = associatedLegendre.calculate(10, 9, -0.25);
        expected = 1.22425639e8;
        break;
    case 9:
        subTestName = std::string("associated_10_11");
        res = associatedLegendre.calculate(10, 11, -1);
        expected = 0;
        break;
    case 10:
        subTestName = std::string("associated_10_-9");
        res = associatedLegendre.calculate(10, -9, 0.3);
        expected = 1.05626898e-9;
        break;
    case 11:
        subTestName = std::string("associated_15_-7");
        res = associatedLegendre.calculate(15, -7, 0.9);
        expected = 1.25421576e-9;
        break;
    case 12:
        subTestName = std::string("associated_100_50");
        res = associatedLegendre.calculate(100, 50, -1.0);
        expected = 0;
        break;
    case 13:
        subTestName = std::string("associated_100_-80");
        res = associatedLegendre.calculate(100, -80, 0.5);
        expected = -4.23304397e-157;
        break;
    case 14:
        subTestName = std::string("associated_1000_100");
        res = associatedLegendre.calculate(1000, 100, -0.25);
        expected = 2.2473217e298;
        break;
    case 15:
        subTestName = std::string("associated_10000_-50");
        res = associatedLegendre.calculate(10000, -50, -0.1);
        expected = 7.16389192e-203;
        break;
    default:
        check(false, "");
        break;
    }
}

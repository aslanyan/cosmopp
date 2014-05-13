#include <test_unit_conversions.hpp>
#include <unit_conversions.hpp>

std::string
TestUnitConversions::name() const
{
    return std::string("UNIT CONVERSIONS");
}

unsigned int
TestUnitConversions::numberOfSubtests() const
{
    return 11;
}

void
TestUnitConversions::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 11, "invalid index " << i);

    using namespace Phys;
    
    switch(i)
    {
    case 0:
        subTestName = std::string("MpcToM");
        res = MpcToM(3.5);
        expected = 1.0799845e23;
        break;
    case 1:
        subTestName = std::string("evToKg");
        res = evToKg(5.0);
        expected = 8.91331e-36;
        break;
    case 2:
        subTestName = std::string("inverseSecToEv");
        res = inverseSecToEv(1.0);
        expected = 6.58212e-16;
        break;
    case 3:
        subTestName = std::string("secToYear");
        res = secToYear(1000);
        expected = 3.16888e-05;
        break;
    case 4:
        subTestName = std::string("unitlessToSec");
        res = unitlessToSec(1.1e20);
        expected = 2.97304e-23;
        break;
    case 5:
        subTestName = std::string("secToUnitless");
        res = secToUnitless(25);
        expected = 9.24978e+43;
        break;
    case 6:
        subTestName = std::string("unitlessToInverseSec");
        res = unitlessToInverseSec(3.3);
        expected = 1.22097e+43;
        break;
    case 7:
        subTestName = std::string("inverseSecToUnitless");
        res = inverseSecToUnitless(15);
        expected = 4.05415e-42;
        break;
    case 8:
        subTestName = std::string("mToUnitless");
        res = mToUnitless(1e-10);
        expected = 1.23416e+24;
        break;
    case 9:
        subTestName = std::string("unitlessToM");
        res = unitlessToM(100);
        expected = 8.10269e-33;
        break;
    case 10:
        subTestName = std::string("kelvinToUnitless");
        res = kelvinToUnitless(2.7);
        expected = 9.55388e-32;
        break;
    default:
        check(false, "");
        break;
    }
}


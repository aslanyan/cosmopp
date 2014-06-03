#include <macros.hpp>
#include <test_int_operations.hpp>
#include <int_operations.hpp>

std::string
TestIntOperations::name() const
{
    return std::string("INT OPERATIONS");
}

unsigned int
TestIntOperations::numberOfSubtests() const
{
    return 2;
}

void
TestIntOperations::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 2, "invalid index " << i);
    
    using namespace Math;
    switch(i)
    {
    case 0:
        subTestName = std::string("isqrt_whole");
        res = isqrt(625);
        expected = 25;
        break;
    case 1:
        subTestName = std::string("isqrt_partial");
        res = isqrt(277829);
        expected = 527;
        break;
    default:
        check(false, "");
        break;
    }
}

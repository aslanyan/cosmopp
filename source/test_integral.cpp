#include <test_integral.hpp>
#include <integral.hpp>

class TestQuadraticFunction : public Math::RealFunction
{
public:
    ~TestQuadraticFunction() {}

    double evaluate(double x) const
    {
        return x * x;
    }
};

std::string
TestIntegral::name() const
{
    return std::string("INTEGRAL");
}

unsigned int
TestIntegral::numberOfSubtests() const
{
    return 1;
}

void
TestIntegral::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);
    
    using namespace Math;

    TestQuadraticFunction f;

    switch(i)
    {
    case 0:
        subTestName = std::string("quadratic");
        res = realIntegral1D(f, -1.0, 3.0, 1000);
        expected = double(28) / 3;
        break;
    default:
        check(false, "");
        break;
    }
}

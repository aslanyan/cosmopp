#include <test_conjugate_gradient.hpp>
#include <conjugate_gradient.hpp>

std::string
TestConjugateGradient::name() const
{
    return std::string("CONJUGATE GRADIENT TESTER");
}

unsigned int
TestConjugateGradient::numberOfSubtests() const
{
    return 1;
}

void
TestConjugateGradient::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);
    
    using namespace Math;
    BasicCGTreats t1(2);
    std::vector<double> b1(2);
    b1[0] = 3.0;
    b1[1] = 4.0;
    ConjugateGradient<BasicCGTreats> cg1(2, &t1, b1);
    switch(i)
    {
    case 0:
        t1.setMatrix(0, 1, 2.0);
        t1.setMatrix(1, 0, 2.0);
        subTestName = std::string("2x2");
        res = cg1.solve()[1];
        expected = double(2) / 3;
        break;
    default:
        check(false, "");
        break;
    }
}

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
    return 4;
}

void
TestConjugateGradient::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < numberOfSubtests(), "invalid index " << i);
    
    using namespace Math;
    BasicCGTreats t1(2);
    std::vector<double> b1(2);
    b1[0] = 3.0;
    b1[1] = 4.0;
    ConjugateGradient<BasicCGTreats> cg1(2, &t1, b1);

    int n2 = 100;
    BasicCGTreats t2(n2);
    std::vector<double> x2(n2), b2(n2);
    for(int i = 0; i < n2; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            const double val = double(n2) - (i - j);
            t2.setMatrix(i, j, val);
            t2.setMatrix(j, i, val);
        }
        x2[i] = i / 5;
    }

    t2.multiplyByMatrix(x2, b2);
    ConjugateGradient<BasicCGTreats> cg2(n2, &t2, b2);

    int n3 = 1000;
    BasicCGTreats t3(n3);
    std::vector<double> x3(n3), b3(n3);
    for(int i = 0; i < n3; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            const double val = double(n3) - (i - j);
            t3.setMatrix(i, j, val);
            t3.setMatrix(j, i, val);
        }
        x3[i] = - double(i) / 10;
    }

    t3.multiplyByMatrix(x3, b3);
    ConjugateGradient<BasicCGTreats> cg3(n3, &t3, b3);

    int n4 = 2500;
    BasicCGTreats t4(n4);
    std::vector<double> x4(n4), b4(n4);
    for(int i = 0; i < n4; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            const double val = double(n4) - (i - j);
            t4.setMatrix(i, j, val);
            t4.setMatrix(j, i, val);
        }
        x4[i] = - double(i) / 10;
    }

    t4.multiplyByMatrix(x4, b4);
    ConjugateGradient<BasicCGTreats> cg4(n4, &t4, b4);
    switch(i)
    {
    case 0:
        t1.setMatrix(0, 1, 2.0);
        t1.setMatrix(1, 0, 2.0);
        subTestName = std::string("2x2");
        res = cg1.solve()[1];
        expected = double(2) / 3;
        break;
    case 1:
        res = cg2.solve()[n2 / 2];
        expected = x2[n2 / 2];
        subTestName = std::string("100x100");
        break;
    case 2:
        res = cg3.solve()[n3 / 2];
        expected = x3[n3 / 2];
        subTestName = std::string("1000x1000");
        break;
    case 3:
        res = cg4.solve()[n4 - 1];
        expected = x4[n4 - 1];
        subTestName = std::string("2500x2500");
        break;
    default:
        check(false, "");
        break;
    }
}

#include <test_polynomial.hpp>
#include <polynomial.hpp>

std::string
TestPolynomial::name() const
{
    return std::string("POLYNOMIAL TESTER");
}

unsigned int
TestPolynomial::numberOfSubtests() const
{
    return 8;
}

void
TestPolynomial::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 8, "invalid index " << i);
    
    using namespace Math;
    Polynomial p(3), q(4);
    p.parameter(0) = 1;
    p.parameter(1) = 2;
    p.parameter(2) = 5;

    q.parameter(0) = 2;
    q.parameter(1) = 7;
    q.parameter(2) = 4;
    q.parameter(3) = -2;
    Polynomial r(1);
    switch(i)
    {
    case 0:
        subTestName = std::string("evaluate");
        res = p.evaluate(3);
        expected = 1 + 2 * 3 + 5 * 3 * 3;
        break;
    case 1:
        subTestName = std::string("multiply_by_number_left");
        r = -3 * p;
        res = r.parameter(1);
        expected = -3 * p.parameter(1);
        break;
    case 2:
        subTestName = std::string("multiply_by_number_right");
        r = p * 2.5;
        res = r.parameter(2);
        expected = 2.5 * p.parameter(2);
        break;
    case 3:
        subTestName = std::string("divide_by_number");
        r = p / 4;
        res = r.parameter(0);
        expected = p.parameter(0) / 4;
        break;
    case 4:
        subTestName = std::string("add");
        r = p + q;
        res = r.parameter(2);
        expected = p.parameter(2) + q.parameter(2);
        break;
    case 5:
        subTestName = std::string("subtract");
        r = p - q;
        res = r.parameter(3);
        expected = -q.parameter(3);
        break;
    case 6:
        subTestName = std::string("multiply");
        r = p * q;
        res = r.parameter(4);
        expected = 16;
        break;
    case 7:
        subTestName = std::string("complex");
        r = (3 * q - p / 2) * p;
        res = r.parameter(2);
        expected = 5 * 5.5 + 20 * 2 + 9.5;
        break;
    default:
        check(false, "");
        break;
    }
}

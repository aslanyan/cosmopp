#include <macros.hpp>
#include <expressions.hpp>

using namespace Expressions;

namespace
{

bool testJustScalar()
{

    JustScalar *xPtr = new JustScalar();
    ScalarPtr x(xPtr);
    const double c = 10.0;
    xPtr->set(c);
    const double v = x->value();
    return v == c;
}

bool testJustScalarDeriv()
{

    JustScalar *xPtr = new JustScalar();
    ScalarPtr x(xPtr);
    JustScalar *yPtr = new JustScalar();
    ScalarPtr y(yPtr);

    const double c1 = 10.0, c2 = -0.5;
    xPtr->set(c1);
    yPtr->set(c2);

    auto d1 = x->derivative(x);
    auto d2 = x->derivative(y);

    if(d1->value() != 1.0)
        return false;
    if(d2->value() != 0.0)
        return false;

    return true;
}

bool testScalarAddSubtract()
{
    JustScalar *xPtr = new JustScalar();
    ScalarPtr x(xPtr);
    JustScalar *yPtr = new JustScalar();
    ScalarPtr y(yPtr);

    auto z = x + y;
    //auto z = ScalarPtr(new ScalarAdd(x, y));

    xPtr->set(10);
    yPtr->set(15);
    
    if(z->value() != 25)
        return false;

    auto zx(z->derivative(x));
    if(zx->value() != 1)
        return false;

    auto diff = x - y;
    if(diff->value() != -5)
        return false;

    auto diffY(diff->derivative(y));
    if(diffY->value() != -1)
        return false;

    auto diffX(diff->derivative(x));
    if(diffX->value() != 1)
        return false;

    return true;
}

bool testScalarOperations()
{
    JustScalar *xPtr = new JustScalar();
    ScalarPtr x(xPtr);
    JustScalar *yPtr = new JustScalar();
    ScalarPtr y(yPtr);

    auto u = x * x + y;
    auto v = 2 * pow(y, 3) + pow(x, 2);
    auto z = u * exp(v);

    const double xVal = 2, yVal = 3;
    const double uVal = xVal * xVal + yVal;
    const double vVal = 2 * std::pow(yVal, 3) + std::pow(xVal, 2);
    const double zVal = uVal * std::exp(vVal);

    xPtr->set(xVal);
    yPtr->set(yVal);

    if(z->value() != zVal)
    {
        output_screen("Value failed!" << std::endl);
        return false;
    }

    auto zx = z->derivative(x);
    const double zxVal = 2 * xVal * std::exp(vVal) * (1 + uVal);
    if(zx->value() != zxVal)
    {
        output_screen("Derivative failed!" << std::endl);
        return false;
    }

    return true;
}

} // namespace

int main()
{
    int res = 0;
    output_screen("Testing Just Scalar..." << std::endl);
    if(testJustScalar())
    {
        ++res;
        output_screen("OK" << std::endl);
    }
    else
    {
        output_screen("FAIL" << std::endl);
    }

    output_screen("Testing Just Scalar Derivative..." << std::endl);
    if(testJustScalarDeriv())
    {
        ++res;
        output_screen("OK" << std::endl);
    }
    else
    {
        output_screen("FAIL" << std::endl);
    }

    output_screen("Testing Scalar Add and Subtract..." << std::endl);
    if(testScalarAddSubtract())
    {
        ++res;
        output_screen("OK" << std::endl);
    }
    else
    {
        output_screen("FAIL" << std::endl);
    }

    output_screen("Testing Scalar Operations Combined..." << std::endl);
    if(testScalarOperations())
    {
        ++res;
        output_screen("OK" << std::endl);
    }
    else
    {
        output_screen("FAIL" << std::endl);
    }

    return res;
}

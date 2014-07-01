#include <table_function.hpp>
#include <test_table_function.hpp>

std::string
TestTableFunction::name() const
{
    return std::string("TABLE FUNCTION TESTER");
}

unsigned int
TestTableFunction::numberOfSubtests() const
{
    return 3;
}

void
TestTableFunction::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 3, "invalid index " << i);

    Math::TableFunction<double, double> t1;
    const double x[3] = {-1, 0, 5};
    const double y[3] = {10, -2, -3};

    for(int i = 0; i < 3; ++i)
        t1[x[i]] = y[i];

    double p;
    switch(i)
    {
    case 0:
        subTestName = std::string("t1");
        p = (x[0] + x[1]) / 2;
        res = t1.evaluate(p);
        expected = (y[0] + y[1]) / 2;
        break;
    case 1:
        subTestName = std::string("t2");
        p = x[1];
        res = t1.evaluate(p);
        expected = y[1];
        break;
    case 2:
        subTestName = std::string("t3");
        p = (x[2] + x[1]) / 2;
        res = t1.evaluate(p);
        expected = (y[2] + y[1]) / 2;
        break;
    default:
        check(false, "");
        break;
    }
}

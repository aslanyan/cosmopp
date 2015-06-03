#include <macros.hpp>
#include <modecode.hpp>
#include <test_modecode.hpp>

std::string
TestModeCode::name() const
{
    return std::string("MODECODE TESTER");
}

unsigned int
TestModeCode::numberOfSubtests() const
{
    return 2;
}

void
TestModeCode::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < numberOfSubtests(), "invalid index " << i);

    subTestName = "quadratic_scalar";
    if(i == 1)
        subTestName = "quadratic_tensor";

    ModeCode::initialize(1, 0.05, 55, true, false, true);
    std::vector<double> vParams(1, -10);

    const bool r = ModeCode::calculate(vParams);
    if(!r)
    {
        output_screen("MODECODE Failed to calculate the power spectra!" << std::endl);
        res = 0;
        expected = 1;
        return;
    }

    res = (i == 0 ? ModeCode::getScalarPs().evaluate(0.05) : ModeCode::getTensorPs().evaluate(0.05));

    // obtained by running modecode by itself.
    expected = (i == 0 ? 5.70997e-09 :7.77754e-10); 
}

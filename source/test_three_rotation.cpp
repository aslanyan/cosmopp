#include <cmath>

#include <three_rotation.hpp>
#include <math_constants.hpp>
#include <test_three_rotation.hpp>

std::string
TestThreeRotation::name() const
{
    return std::string("THREE ROTATION TESTER");
}

unsigned int
TestThreeRotation::numberOfSubtests() const
{
    return 2;
}

void
TestThreeRotation::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < numberOfSubtests(), "invalid index " << i);
    
    using namespace Math;

    ThreeRotationMatrix m1(pi / 4, 0, 0);
    ThreeVectorDouble v1(1, 0, 0);
    ThreeVectorDouble v2(0, 1, 0);
    ThreeVectorDouble r;

    ThreeRotationMatrix m2(0, -pi / 2, 0);

    switch(i)
    {
    case 0:
        subTestName = std::string("rot_z");
        r = m1 * v1;
        res = r.x();
        expected = std::sqrt(2.0) / 2.0;
        break;
    case 1:
        subTestName = std::string("rot_x");
        r = m2 * v2;
        res = r.z();
        expected = 1;
        break;
    default:
        check(false, "");
        break;
    }
}

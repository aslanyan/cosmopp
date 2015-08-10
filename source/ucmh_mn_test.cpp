#include <cmath>

#include <macros.hpp>
#include <mn_scanner.hpp>

namespace
{

class MyLike : public Math::LikelihoodFunction
{
public:
    MyLike() {}
    virtual ~MyLike() {}

    virtual double calculate(double* params, int nParams)
    {
        check(nParams == 2, "");
        check(params[0] != 0, "");
        //const double x = std::pow(10.0, params[1]);
        const double x = params[1];
        const double r = params[0] / x;

        const double rMean = 1.0;
        const double rSigma = 0.005;

        const double diff = (r - rMean);
        return diff * diff / (rSigma * rSigma);
    }
};

}

int main(int argc, char *argv[])
{
    MyLike like;
    std::string root = "mn_ucmh_test";
    MnScanner mn(2, like, 500, root);
    mn.setParam(0, "e", 0, 0.1);
    //mn.setParam(1, "v", -10, 10);
    mn.setParam(1, "v", 0, 0.1);
    mn.run(true);
    return 0;
}

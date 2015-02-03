#include <ctime>
#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <gaussian_process.hpp>
#include <random.hpp>
#include <markov_chain.hpp>
#include <test_gaussian_process.hpp>

std::string
TestGaussianProcess::name() const
{
    return std::string("GAUSSIAN PROCESS TESTER");
}

unsigned int
TestGaussianProcess::numberOfSubtests() const
{
    return 1;
    //return 2;
}

namespace
{

double f(double x)
{
    return 5 * x * x - 3 * x + 10;
}

double f2(double x, double y, double z)
{
    const double x1 = 100 * x + y;
    const double y1 = -2 * y + 10 * z * z;
    const double z1 = z + std::cos(x) - 10 * y;

    return (std::sin(x1) + z1) * y1 * y1 + 3 * y1 * y1;
}

}

void
TestGaussianProcess::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    switch(i)
    {
    case 0:
        runSubTest0(res, expected, subTestName);
        break;
    case 1:
        runSubTest1(res, expected, subTestName);
        break;
    default:
        check(false, "invalid subtest " << i);
        break;
    }
}

void
TestGaussianProcess::runSubTest0(double& res, double& expected, std::string& subTestName)
{
    std::vector<std::vector<double> > points;
    std::vector<double> data;

    const double xMin = -1.0;
    const double xMax = 1.0;

    Math::UniformRealGenerator gen(std::time(0), xMin, xMax);

    const int nPoints = 10;

    std::vector<double> p(1);
    double d;

    for(int i = 0; i < nPoints; ++i)
    {
        p[0] = gen.generate();
        points.push_back(p);
    }

    for(int i = 0; i < points.size(); ++i)
    {
        d = f(points[i][0]);
        data.push_back(d);
    }

    Math::GaussianProcess gp(1);

    gp.set(points, data, true);

    std::vector<std::vector<double> > input(1);
    input[0].resize(1, 0);
    std::vector<double> mean;
    LaGenMatDouble covariance;

    gp.calculate(input, mean, covariance);

    check(covariance(0, 0) >= 0, "");
    const double error = std::sqrt(covariance(0, 0));

    output_screen("Result = " << mean[0] << " +/- " << std::sqrt(covariance(0, 0)) << std::endl);

    subTestName = "simple";

    expected = 1;
    res = 1;
    if(error > 100.0)
    {
        res = 0;
        output_screen("FAIL! Error is " << error << ". Too large!" << std::endl);
    }

    if(mean[0] < f(input[0][0]) - 5 * error || mean[0] > f(input[0][0]) + 5 * error)
    {
        res = 0;
        output_screen("FAIL! The mean value " << mean[0] << " is more than 5 sigma away from the expected value " << f(input[0][0]) << ". (sigma = " << error << ")." << std::endl);
    }

    // writing into files for plotting
    std::ofstream out("test_files/gaussian_process_test_points.txt");
    StandardException exc;
    if(!out)
    {
        std::string exceptionStr = "Cannot write into file test_files/gaussian_process_test_points.txt!";
        exc.set(exceptionStr);
        throw exc;
    }

    for(int i = 0; i < nPoints; ++i)
        out << points[i][0] << '\t' << data[i] << std::endl;

    out.close();
    out.open("test_files/gaussian_process_test_results.txt");
    if(!out)
    {
        std::string exceptionStr = "Cannot write into file test_files/gaussian_process_test_results.txt!";
        exc.set(exceptionStr);
        throw exc;
    }

    const int nTest = 1000;
    const double delta = (xMax - xMin) / nTest;
    input.clear();
    input.resize(nTest + 1);
    for(int i = 0; i <= nTest; ++i)
    {
        const double x = xMin + i * delta;
        input[i].resize(1, x);
    }
    gp.calculate(input, mean, covariance);

    for(int i = 0; i <= nTest; ++i)
    {
        const double e = std::sqrt(covariance(i, i));
        out << input[i][0] << '\t' << f(input[i][0]) << '\t' << mean[i] << '\t' << e << '\t' << mean[i] - e << '\t' << mean[i] + e << std::endl;
    }
    out.close();
}

void
TestGaussianProcess::runSubTest1(double& res, double& expected, std::string& subTestName)
{
    subTestName = "2d";
    res = 1;
    expected = 1;

    const int n = 2000;

    std::vector<std::vector<double> > points, testPoints;
    std::vector<double> data, testData;

    std::vector<double> p(3);
    double d;

    Math::UniformRealGenerator gen(std::time(0), -1, 1);

    for(int i = 0; i < n; ++i)
    {
        p[0] = gen.generate();
        p[1] = gen.generate();
        p[2] = gen.generate();

        d = f2(p[0], p[1], p[2]);

        if(i % 2)
        {
            points.push_back(p);
            data.push_back(d);
        }
        else
        {
            testPoints.push_back(p);
            testData.push_back(d);
        }
    }

    Math::GaussianProcess gp(3);
    gp.set(points, data, true);
    std::vector<double> mean;
    LaGenMatDouble covariance;

    Posterior1D post;
    for(int i = 0; i < testPoints.size(); ++i)
    {
        std::vector<std::vector<double> > in;
        in.push_back(testPoints[i]);
        std::vector<double> mean;
        LaGenMatDouble covariance;

        gp.calculate(in, mean, covariance);

        const double estimatedError = std::sqrt(covariance(0, 0));

        const double correctError = std::abs(mean[0] - testData[i]);

        if(estimatedError == 0)
        {
            check(correctError == 0, "");
        }
        else
        {
            const double ratio = correctError / estimatedError;
            post.addPoint(ratio, 1, 1);
        }

        //output_screen1("Result = " << mean[0] << " +/- " << estimatedError << ", correct result = " << testData[i] << std::endl);
    }

    post.generate();
    output_screen1("Posterior 1 sigma is: " << post.get1SigmaUpper() << std::endl);
    post.writeIntoFile("test_files/gaussian_process_test_error_ratio.txt");
}

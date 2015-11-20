#include <cosmo_mpi.hpp>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>

#include <pde.hpp>
#include <macros.hpp>
#include <math_constants.hpp>
#include <numerics.hpp>
#include <timer.hpp>
#include <test_pde.hpp>

std::string
TestPDE::name() const
{
    return std::string("PDE TESTER");
}

unsigned int
TestPDE::numberOfSubtests() const
{
    return 6;
}

namespace
{

class PDE0 : public Math::InitialValPDEInterface
{
public:
    PDE0(bool source = false): source_(source) {}
    virtual ~PDE0() {}

    virtual int spaceDim() const { return 1; }
    virtual int funcDim() const { return 1; }

    virtual void evaluate(double t, const double *x, const std::vector<double>& u, std::vector<std::vector<double> > *f, std::vector<double> *s) const
    {
        check(u.size() == 1, "");
        check(f->size() == 1, "");
        check((*f)[0].size() == 1, "");
        check(s->size() == 1, "");
        (*f)[0][0] = u[0];
        (*s)[0] = (source_ ? ((x[0] - 2 * Math::pi) * (x[0] - 2 * Math::pi) - 2 * t * (x[0] - 2 * Math::pi)) / 100. : 0.0);
    }
private:
    bool source_;
};

class PDE2 : public Math::InitialValPDEInterface
{
public:
    PDE2() {}
    virtual ~PDE2() {}

    virtual int spaceDim() const { return 1; }
    virtual int funcDim() const { return 3; }

    virtual void evaluate(double t, const double *x, const std::vector<double>& u, std::vector<std::vector<double> > *f, std::vector<double> *s) const
    {
        check(u.size() == 3, "");
        check(f->size() == 1, "");
        check((*f)[0].size() == 3, "");
        check(s->size() == 3, "");
        (*f)[0][0] = 0;
        (*f)[0][1] = u[2];
        (*f)[0][2] = u[1];
        (*s)[0] = u[2];
        (*s)[1] = 0;
        (*s)[2] = 0;
    }
};

class PDE3 : public Math::InitialValPDEInterface
{
public:
    PDE3() {}
    virtual ~PDE3() {}

    virtual int spaceDim() const { return 1; }
    virtual int funcDim() const { return 2; }

    virtual void evaluate(double t, const double *x, const std::vector<double>& u, std::vector<std::vector<double> > *f, std::vector<double> *s) const
    {
        check(u.size() == 2, "");
        check(f->size() == 1, "");
        check((*f)[0].size() == 2, "");
        check(s->size() == 2, "");
        (*f)[0][0] = u[0];
        (*f)[0][1] = u[1];
        (*s)[0] = -2. * u[1];
        (*s)[1] = 2. * u[0];
    }
};

class PDE4 : public Math::InitialValPDEInterface
{
public:
    PDE4() {}
    virtual ~PDE4() {}

    virtual int spaceDim() const { return 2; }
    virtual int funcDim() const { return 1; }

    virtual void evaluate(double t, const double *x, const std::vector<double>& u, std::vector<std::vector<double> > *f, std::vector<double> *s) const
    {
        check(u.size() == 1, "");
        check(f->size() == 2, "");
        check((*f)[0].size() == 1, "");
        check((*f)[1].size() == 1, "");
        check(s->size() == 1, "");
        (*f)[0][0] = 0;
        (*f)[1][0] = 0;
        (*s)[0] = -std::sin(t + x[0] + x[1]);
    }
};

class PDE5 : public Math::InitialValPDEInterface
{
public:
    PDE5() {}
    virtual ~PDE5() {}

    virtual int spaceDim() const { return 2; }
    virtual int funcDim() const { return 1; }

    virtual void evaluate(double t, const double *x, const std::vector<double>& u, std::vector<std::vector<double> > *f, std::vector<double> *s) const
    {
        check(u.size() == 1, "");
        check(f->size() == 2, "");
        check((*f)[0].size() == 1, "");
        check((*f)[1].size() == 1, "");
        check(s->size() == 1, "");
        (*f)[0][0] = u[0];
        (*f)[1][0] = u[0];
        (*s)[0] = std::sin(t - x[0] + 3. * x[1]);
    }
};

class PDE0InitialFunc : public Math::RealFunctionMultiDim
{
public:
    PDE0InitialFunc() {}
    virtual ~PDE0InitialFunc() {}
    virtual double evaluate(const std::vector<double>& x) const { check(x.size() == 1, ""); return std::cos(x[0]); }
};

class PDE2InitialFunc0 : public Math::RealFunctionMultiDim
{
public:
    PDE2InitialFunc0() {}
    virtual ~PDE2InitialFunc0() {}
    virtual double evaluate(const std::vector<double>& x) const { check(x.size() == 1, ""); return std::cos(x[0]); }
};

class PDE2InitialFunc1 : public Math::RealFunctionMultiDim
{
public:
    PDE2InitialFunc1() {}
    virtual ~PDE2InitialFunc1() {}
    virtual double evaluate(const std::vector<double>& x) const { check(x.size() == 1, ""); return -std::sin(x[0]); }
};

class PDE3InitialFunc0 : public Math::RealFunctionMultiDim
{
public:
    PDE3InitialFunc0() {}
    virtual ~PDE3InitialFunc0() {}
    virtual double evaluate(const std::vector<double>& x) const { check(x.size() == 1, ""); return std::cos(-x[0]); }
};

class PDE3InitialFunc1 : public Math::RealFunctionMultiDim
{
public:
    PDE3InitialFunc1() {}
    virtual ~PDE3InitialFunc1() {}
    virtual double evaluate(const std::vector<double>& x) const { check(x.size() == 1, ""); return std::sin(-x[0]); }
};

class PDE4InitialFunc : public Math::RealFunctionMultiDim
{
public:
    PDE4InitialFunc() {}
    virtual ~PDE4InitialFunc() {}
    virtual double evaluate(const std::vector<double>& x) const { check(x.size() == 2, ""); return std::cos(x[0] + x[1]); }
};

class PDE5InitialFunc : public Math::RealFunctionMultiDim
{
public:
    PDE5InitialFunc() {}
    virtual ~PDE5InitialFunc() {}
    virtual double evaluate(const std::vector<double>& x) const { check(x.size() == 2, ""); return std::cos(-x[0] + 3. * x[1]); }
};

}

void
TestPDE::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    if(i < 4)
    {
        res = 1;
        expected = 1;
        subTestName = "not_run";
        //return;
    }

    std::unique_ptr<Math::InitialValPDEInterface> pde;

    // set subTestName and pde
    switch(i)
    {
    case 0:
        subTestName = "simplest_1d";
        pde.reset(new PDE0(false));
        break;
    case 1:
        subTestName = "simplest_1d_with_source";
        pde.reset(new PDE0(true));
        break;
    case 2:
        subTestName = "not_so_simple_1d";
        pde.reset(new PDE2);
        break;
    case 3:
        subTestName = "another_1d";
        pde.reset(new PDE3);
        break;
    case 4:
        subTestName = "simplest_2d";
        pde.reset(new PDE4);
        break;
    case 5:
        subTestName = "simple_2d";
        pde.reset(new PDE5);
        break;
    default:
        check(false, "");
        break;
    }

    Math::InitialValPDESolver solver(&(*pde));

    std::vector<std::unique_ptr<Math::RealFunctionMultiDim> > initialFuncs(pde->funcDim());
    std::vector<Math::RealFunctionMultiDim*> initial(pde->funcDim());

    std::vector<double> xMin(pde->spaceDim()), xMax(pde->spaceDim());
    std::vector<int> nx(pde->spaceDim());
    double tMax;

    // set initial, xMin, xMax, nx, tMax
    switch(i)
    {
    case 0:
    case 1:
        initialFuncs[0].reset(new PDE0InitialFunc);
        xMin[0] = 0;
        xMax[0] = 4 * Math::pi;
        nx[0] = 512;
        tMax = 10;
        break;
    case 2:
        initialFuncs[0].reset(new PDE2InitialFunc0);
        initialFuncs[1].reset(new PDE2InitialFunc1);
        initialFuncs[2].reset(new PDE2InitialFunc1);
        xMin[0] = 0;
        xMax[0] = 4 * Math::pi;
        nx[0] = 512;
        tMax = 10;
        break;
    case 3:
        initialFuncs[0].reset(new PDE3InitialFunc0);
        initialFuncs[1].reset(new PDE3InitialFunc1);
        xMin[0] = -2 * Math::pi;
        xMax[0] = 2 * Math::pi;
        nx[0] = 512;
        tMax = 10;
        break;
    case 4:
        initialFuncs[0].reset(new PDE4InitialFunc);
        xMin[0] = 0;
        xMax[0] = 2 * Math::pi;
        xMin[1] = -2 * Math::pi;
        xMax[1] = 2 * Math::pi;
        nx[0] = 128;
        nx[1] = 256;
        tMax = 10;
        break;
    case 5:
        initialFuncs[0].reset(new PDE5InitialFunc);
        xMin[0] = 0;
        xMax[0] = 2 * Math::pi;
        xMin[1] = 0;
        xMax[1] = 2 * Math::pi;
        nx[0] = 128;
        nx[1] = 128;
        tMax = 10;
        break;
    default:
        check(false, "");
        break;
    }

    for(int i = 0; i < pde->funcDim(); ++i)
        initial[i] = &(*initialFuncs[i]);

    solver.set(initial, xMin, xMax, nx);

    std::stringstream timerName;
    timerName << "PDE" << i;
    Timer t1(timerName.str().c_str());
    t1.start();
    solver.propagate(tMax);
    t1.end();
    const double t = solver.getCurrentT();

    check(Math::areEqual(tMax, t, 1e-2), "");

    res = 1;
    expected = 1;

    int dimToTest = rand_.generate() % pde->spaceDim();
    int fieldToTest = rand_.generate() % pde->funcDim();
    int pointToTest = rand_.generate() % solver.getNx(dimToTest);

    for(int j = 0; j < pde->funcDim(); ++j)
    {
        std::stringstream outputFileName;
        outputFileName << "test_files/pde" << i << '_' << j;
        if(CosmoMPI::create().numProcesses() > 1)
            outputFileName << '_' << CosmoMPI::create().processId();
        outputFileName << ".txt";
        std::ofstream outFile(outputFileName.str().c_str());

        if(pde->spaceDim() == 1)
        {
            std::vector<double> coords(1);
            std::vector<int> ind(1);
            for(ind[0] = -1; ind[0] <= solver.getNx(0); ++ind[0])
            {
                solver.physicalCoords(&(ind[0]), &(coords[0]));
                const double x = coords[0];
                const double y = solver.getField(ind)[j];
                double expY;
                // set expected value
                switch(i)
                {
                case 0:
                    check(j == 0, "");
                    expY = std::cos(t + x);
                    break;
                case 1:
                    check(j == 0, "");
                    expY = std::cos(t + x) + t * (x - 2 * Math::pi) * (x - 2 * Math::pi) / 100.;
                    break;
                case 2:
                    check(j < 3, "");
                    switch(j)
                    {
                    case 0:
                        expY = std::cos(t + x);
                        break;
                    case 1:
                    case 2:
                        expY = -std::sin(t + x);
                        break;
                    }
                    break;
                case 3:
                    check(j < 2, "");
                    expY = (j == 0 ? std::cos(t - x) : std::sin(t - x));
                    break;
                default:
                    check(false, "");
                    break;
                }
                outFile << x << '\t' << y << '\t' << expY << std::endl;

                check(dimToTest == 0, "");
                if(isMaster() && j == fieldToTest && ind[0] == pointToTest)
                {
                    res = y;
                    expected = expY;
                }
            }
        }
        else if(pde->spaceDim() == 2)
        {
            std::vector<double> coords(2);
            std::vector<int> ind(2);
            for(ind[0] = -1; ind[0] <= solver.getNx(0); ++ind[0])
            {
                for(ind[1] = -1; ind[1] <= solver.getNx(1); ++ind[1])
                {
                    solver.physicalCoords(&(ind[0]), &(coords[0]));
                    const double x = coords[0];
                    const double y = coords[1];
                    const double z = solver.getField(ind)[j];
                    double expZ;
                    // set expected value
                    switch(i)
                    {
                    case 4:
                        check(j == 0, "");
                        expZ = std::cos(t + x + y);
                        break;
                    case 5:
                        check(j == 0, "");
                        expZ = std::cos(t - x + 3. * y);
                        break;
                    default:
                        check(false, "");
                        break;
                    }
                    outFile << x << '\t' << y << '\t' << z << '\t' << expZ << std::endl;
                    if(isMaster() && j == fieldToTest && ind[dimToTest] == pointToTest)
                    {
                        res = z;
                        expected = expZ;
                    }
                }
            }
        }
        else
        {
            check(false, "not implemented");
        }
        outFile.close();
    }
}


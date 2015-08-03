#ifndef COSMO_PP_TAYLOR_PK_HPP
#define COSMO_PP_TAYLOR_PK_HPP

#include <vector>

#include <macros.hpp>
#include <table_function.hpp>

struct precision;
struct background;
struct thermo;
struct perturbs;
struct transfers;
struct primordial;
struct spectra;
struct nonlinear;
struct lensing;
struct output;

class TaylorPk
{
public:
    TaylorPk(double kPivot = 0.05, double kMin = 1e-6, double kMax = 1.0, int kPerDecade = 10);
    ~TaylorPk();

    void addKValue(double k, double sMin = 0, double sMax = 1, double tMin = 0, double tMax = 1);

    bool calculate(const std::vector<double>& vParams, double *badLike = NULL);

    inline const Math::TableFunction<double, double>& getScalarPs() { return scalarPs_; }
    inline const Math::TableFunction<double, double>& getTensorPs() { return tensorPs_; }

private:
    Math::TableFunction<double, double> scalarPs_;
    Math::TableFunction<double, double> tensorPs_;

    std::map<double, double> scalarLower_, scalarUpper_, tensorLower_, tensorUpper_;

    precision* pr_;
    background* br_;
    thermo* th_;
    perturbs* pt_;
    transfers* tr_;
    primordial* pm_;
    spectra* sp_;
    nonlinear* nl_;
    lensing* le_;
    output* op_;
};

#endif

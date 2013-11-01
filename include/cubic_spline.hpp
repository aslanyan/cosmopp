#ifndef COSMO_CPP_CUBIC_SPLINE_HPP
#define COSMO_CPP_CUBIC_SPLINE_HPP

#include <vector>
#include <set>

#include <function.hpp>
#include <macros.hpp>

namespace Math
{

/// A class that implements a cubic spline.
class CubicSpline : public RealFunction
{
public:
    /// Constructor.
    /// \param x A vector containing the x coordinates of the points to be used for the spline. At least 2 points needed. No double entries are allowed. The vector does not have to be sorted.
    /// \param y A vector containing the y coordinates of the points to be used for the spline. Has to be of the same size as x.
    inline CubicSpline(const std::vector<double>& x, const std::vector<double>& y);

    /// Destructor.
    inline ~CubicSpline();

    /// A function that evaluates the spline for any given x. The argument x has to be between the smallest and biggest values of x of all of the spline points.
    inline virtual double evaluate(double x) const;
private:
    struct SplineData
    {
        double a, b, c, d;
    };

    struct SplinePoints
    {
        double x, y;
        SplineData* data;

        bool operator < (const SplinePoints& other) const { return x < other.x; }
    };

    typedef std::set<SplinePoints> SplineDataSet;

private:
    SplineDataSet splineDataSet_;
};

CubicSpline::CubicSpline(const std::vector<double>& x, const std::vector<double>& y)
{
    check(x.size() == y.size(), "the x and y vectors need to have the same size");
    check(x.size() >= 2, "at least 2 points needed for the spline");

    SplinePoints sp;
    for(int i = 0; i < x.size(); ++i)
    {
        sp.x = x[i];
        sp.y = y[i];
        sp.data = new SplineData;
        splineDataSet_.insert(sp);
    }

    check(splineDataSet_.size() == x.size(), "duplicate entries detected");

    // algorithm from http://en.wikipedia.org/w/index.php?title=Spline_%28mathematics%29&oldid=288288033#Algorithm_for_computing_natural_cubic_splines
    std::vector<double> a, h;

    for(SplineDataSet::const_iterator it = splineDataSet_.begin(); it != splineDataSet_.end(); ++it)
        a.push_back((*it).y);

    check(a.size() == x.size(), "");

    SplineDataSet::const_iterator it = splineDataSet_.begin();
    double xPrev = (*it).x;
    while(++it != splineDataSet_.end())
    {
        double x = (*it).x;
        h.push_back(x - xPrev);
        xPrev = x;
    }

    check(h.size() == a.size() - 1, "");

    std::vector<double> alpha(h.size());
    for(int i = 1; i < alpha.size(); ++i)
        alpha[i] = (3.0 / h[i]) * (a[i + 1] - a[i]) - (3.0 / h[i - 1]) * (a[i] - a[i - 1]);

    std::vector<double> c(a.size()), l(a.size()), mu(a.size()), z(a.size());
    std::vector<double> b(h.size()), d(h.size());
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for(int i = 1; i < h.size(); ++i)
    {
        l[i] = 2 * (h[i] + h[i - 1]) - h[i - 1] * mu[i - 1];
        
        check(l[i] != 0, "bad things happening");

        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[l.size() - 1] = 1;
    z[z.size() - 1] = 0;
    c[c.size() - 1] = 0;

    for(int j = h.size() - 1; j >= 0; --j)
    {
        c[j] = z[j] - mu[j] * c[j + 1];

        check(h[j] != 0, "bad things happening");

        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    SplineDataSet::iterator iter = splineDataSet_.begin();
    for(int i = 0; i < h.size(); ++i)
    {
        check((*iter).y == a[i], "");
        (*iter).data->a = a[i];
        (*iter).data->b = b[i];
        (*iter).data->c = c[i];
        (*iter).data->d = d[i];

        ++iter;
    }

    check(++iter == splineDataSet_.end(), "");
}

CubicSpline::~CubicSpline()
{
    for(SplineDataSet::iterator it = splineDataSet_.begin(); it != splineDataSet_.end(); ++it)
        delete (*it).data;
}

double
CubicSpline::evaluate(double x) const
{
    check(splineDataSet_.size() >= 2, "not properly initialized");
    SplinePoints xx;
    xx.x = x;
    SplineDataSet::const_iterator it = splineDataSet_.lower_bound(xx);
    check(it != splineDataSet_.end(), "element is outside the range");
    if(it == splineDataSet_.begin())
    {
        check((*it).x == x, "element is outside the range");
        return (*it).y;
    }

    --it;
    check(x > (*it).x, "");
    const double deltaX = x - (*it).x;
    const double deltaX2 = deltaX * deltaX;
    const double deltaX3 = deltaX2 * deltaX;

    const SplineData& sp = *((*it).data);
    return sp.a + sp.b * deltaX + sp.c * deltaX2 + sp.d * deltaX3;
}

} // namespace Math

#endif


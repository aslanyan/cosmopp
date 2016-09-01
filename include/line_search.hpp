#ifndef COSMO_PP_LINE_SEARCH_HPP
#define COSMO_PP_LINE_SEARCH_HPP

#include <algorithm>
#include <cmath>

#include <macros.hpp>

namespace Math
{

inline
int moreThuenteStep(double &stx, double &fx, double &dx, double &sty, double &fy, double &dy, double &stp, double &fp, double& dp, bool &brackt, double stpmin, double stpmax)
{
    const double p66 = 0.66;
    int info = 0;

    if((brackt && (stp <= std::min(stx, sty) || stp >= std::max(stx, sty))) || dx * (stp - stx) >= 0 || stpmax < stpmin)
        return info;
    check(dx != 0, "");
    double sgnd = dp * (dx / std::abs(dx));

    bool bound = false;
    double theta, s, gamma, p, q, r, stpc = 0, stpq = 0, stpf = 0;

    if(fp > fx)
    {
        info = 1;
        bound = true;
        theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        s = std::max(std::max(std::abs(theta), std::abs(dx)), std::abs(dp));
        check(s > 0, "");
        gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if(stp < stx)
            gamma = -gamma;
        p = (gamma - dx) + theta;
        q = ((gamma - dx) + gamma) + dp;
        r = p / q;
        stpc = stx + r * (stp - stx);
        stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx);
        if(std::abs(stpc - stx) < std::abs(stpq - stx))
            stpf = stpc;
        else
            stpf = stpc + (stpq - stpc) / 2;
        brackt = true;
    }
    else if(sgnd < 0)
    {
        info = 2;
        bound = false;
        theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        s = std::max(std::max(std::abs(theta), std::abs(dx)), std::abs(dp));
        check(s > 0, "");
        gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if(stp > stx)
            gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dx;
        r = p / q;
        stpc = stp + r * (stx - stp);
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if(std::abs(stpc - stp) > std::abs(stpq - stp))
            stpf = stpc;
        else
            stpf = stpq;
        brackt = true;
    }
    else if(std::abs(dp) < std::abs(dx))
    {
        info = 3;
        bound = true;
        theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        s = std::max(std::max(std::abs(theta), std::abs(dx)), std::abs(dp));
        check(s > 0, "");
        gamma = s * std::sqrt(std::max(0.0, (theta / s) * (theta / s) - (dx / s) * (dp / s)));
        if(stp > stx)
            gamma = -gamma;
        p = (gamma - dp) + theta;
        q = (gamma + (dx - dp)) + gamma;
        r = p / q;
        if(r < 0 && gamma != 0)
            stpc = stp + r * (stx - stp);
        else if(stp > stx)
            stpc = stpmax;
        else
            stpc = stpmin;
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if(brackt)
        {
            if(std::abs(stp - stpc) < std::abs(stp - stpq))
                stpf = stpc;
            else
                stpf = stpq;
        }
        else
        {
            if(std::abs(stp - stpc) > std::abs(stp - stpq))
                stpf = stpc;
            else
                stpf = stpq;
        }
    }
    else
    {
        info = 4;
        bound = false;
        if(brackt)
        {
            theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
            s = std::max(std::max(std::abs(theta), std::abs(dy)), std::abs(dp));
            check(s > 0, "");
            gamma = s * std::sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
            if(stp > sty)
                gamma = -gamma;
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy;
            r = p / q;
            stpc = stp + r * (sty - stp);
            stpf = stpc;
        }
        else if(stp > stx)
            stpf = stpmax;
        else
            stpf = stpmin;
    }

    if(fp > fx)
    {
        sty = stp;
        fy = fp;
        dy = dp;
    }
    else
    {
        if(sgnd < 0)
        {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }

    stpf = std::min(stpmax, stpf);
    stpf = std::max(stpmin, stpf);
    stp = stpf;
    if(brackt && bound)
    {
        if(sty > stx)
            stp = std::min(stx + p66 * (sty - stx), stp);
        else
            stp = std::max(stx + p66 * (sty - stx), stp);
    }
    
    return info;
}

template <typename LargeVector, typename Function>
int moreThuenteSearch(Function *func, const LargeVector& x0, double &f, const LargeVector& g0, const LargeVector& s, double &stp, double ftol, double gtol, double xtol, double stpmin, double stpmax, int maxfev, LargeVector *x, LargeVector *g, int &nfev)
{
    check(stp > 0, "");
    check(ftol >= 0, "");
    check(gtol >= 0, "");
    check(xtol >= 0, "");
    check(stpmin >= 0, "");
    check(stpmax >= stpmin, "");
    check(maxfev > 0, "");

    x->copy(x0);
    g->copy(g0);

    const double dginit = g->dotProduct(s);
    check(dginit < 0, "");

    const double p66 = 0.66;
    int xtrapf = 4;
    int info = 0;
    int infoc = 1;
    bool brackt = false;
    bool stage1 = true;

    nfev = 0;
    double finit = f;
    double dgtest = ftol * dginit;
    double width = stpmax - stpmin;
    double width1 = 2 * width;

    double stx = 0;
    double fx = finit;
    double dgx = dginit;
    double sty = 0;
    double fy = finit;
    double dgy = dginit;

    double stmin, stmax, dg, ftest1, fm, fxm, fym, dgm, dgxm, dgym;
    while(true)
    {
        if(brackt)
        {
            stmin = std::min(stx, sty);
            stmax = std::max(stx, sty);
        }
        else
        {
            stmin = stx;
            stmax = stp + xtrapf * (stp - stx);
        }

        stp = std::max(stp, stpmin);
        stp = std::min(stp, stpmax);

        if((brackt && (stp <= stmin || stp >= stmax)) || nfev >= maxfev - 1 || infoc == 0 || (brackt && stmax - stmin <= xtol * stmax))
            stp = stx;

        x->copy(x0);
        x->add(s, stp);
        func->set(*x);
        f = func->value();
        func->derivative(g);
        ++nfev;
        dg = g->dotProduct(s);
        ftest1 = finit + stp * dgtest;

        if((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0)
            info = 6;
        if(stp == stpmax && f <= ftest1 && dg <= dgtest)
            info = 5;
        if(stp == stpmin && (f > ftest1 || dg >= dgtest))
            info = 4;
        if(nfev >= maxfev)
            info = 3;
        if(brackt && stmax - stmin <= xtol * stmax)
            info = 2;
        if(f <= ftest1 && std::abs(dg) <= gtol * (-dginit))
            info = 1;

        if(info != 0)
            return info;

        if(stage1 && f <= ftest1 && dg >= std::min(ftol, gtol) * dginit)
            stage1 = false;

        if(stage1 && f <= fx && f > ftest1)
        {
            fm = f - stp * dgtest;
            fxm = fx - stx * dgtest;
            fym = fy - sty * dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;

            infoc = moreThuenteStep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax);

            fx = fxm + stx * dgtest;
            fy = fym + sty * dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
        }
        else
            infoc = moreThuenteStep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax);

        if(brackt)
        {
            if(std::abs(sty - stx) >= p66 * width1)
                stp = stx + (sty - stx) / 2;
            width1 = width;
            width = std::abs(sty - stx);
        }
    }
}

} // namespace Math

#endif


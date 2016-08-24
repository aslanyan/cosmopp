#ifndef COSMO_PP_EXPRESSIONS_HPP
#define COSMO_PP_EXPRESSIONS_HPP

#include <memory>
#include <cmath>

#include <macros.hpp>

namespace Expressions
{

class Scalar;
class Vector;
class Vector2Vector;

using ScalarPtr = typename std::shared_ptr<Scalar>;
using VectorPtr = typename std::shared_ptr<Vector>;
using Vector2VectorPtr = typename std::shared_ptr<Vector2Vector>;

class Scalar
{
public:
    virtual ~Scalar() {}

    virtual double value() = 0;

    virtual ScalarPtr derivative(ScalarPtr other) = 0;
    virtual VectorPtr derivative(VectorPtr other) = 0;
};

class Vector
{
public:
    virtual ~Vector() {}

    virtual double sum() = 0;

    virtual VectorPtr derivative(ScalarPtr other) = 0;
    virtual Vector2VectorPtr derivative(VectorPtr other) = 0;
};

class InScalar
{
public:
    virtual ~InScalar() {}

    void setInputScalar(ScalarPtr n) { inN_ = n; }

protected:
    ScalarPtr inN_;
};

class InScalarScalar
{
public:
    virtual ~InScalarScalar() {}

    void setInputScalar1(ScalarPtr n) { inN1_ = n; }
    void setInputScalar2(ScalarPtr n) { inN2_ = n; }

protected:
    ScalarPtr inN1_;
    ScalarPtr inN2_;
};

class InVector
{
public:
    virtual ~InVector() {}

    virtual void setInputVector(VectorPtr v) = 0;
};

class InVectorVector
{
public:
    virtual ~InVectorVector() {}

    virtual void setInputVector1(VectorPtr v) = 0;
    virtual void setInputVector2(VectorPtr v) = 0;
};

class InScalarVector
{
public:
    virtual ~InScalarVector() {}

    virtual void setInputScalar(ScalarPtr n) = 0;
    virtual void setInputVector(VectorPtr v) = 0;
};

class Scalar2Scalar : public Scalar, public InScalar
{
public:
    virtual ~Scalar2Scalar() {}

    // Scalar functions
    virtual double value() = 0;
    virtual ScalarPtr derivative(ScalarPtr other) = 0;
    virtual VectorPtr derivative(VectorPtr other) = 0;
};

class ScalarScalar2Scalar : public Scalar, public InScalarScalar
{
public:
    virtual ~ScalarScalar2Scalar() {}

    // Scalar functions
    virtual double value() = 0;
    virtual ScalarPtr derivative(ScalarPtr other) = 0;
    virtual VectorPtr derivative(VectorPtr other) = 0;
};

class Vector2Scalar : public Scalar, public InVector
{
public:
    virtual ~Vector2Scalar() {}
};

class Vector2Vector : public Vector, public InVector
{
public:
    virtual ~Vector2Vector() {}
};

class VectorVector2Vector : public Vector, public InVectorVector
{
public:
    virtual ~VectorVector2Vector() {}
};

class ZeroVector : public Vector
{
public:
    ZeroVector() {}
    virtual ~ZeroVector() {}

    double sum() { return 0; }

    virtual VectorPtr derivative(ScalarPtr other) { return VectorPtr(new ZeroVector); }
    virtual Vector2VectorPtr derivative(VectorPtr other)
    {
        //TBD
        return Vector2VectorPtr();
    }
};


class JustScalar : public Scalar
{
public:
    JustScalar(double val = 0) : val_(val) {}
    virtual ~JustScalar() {}

    void set(double val) { val_ = val; }
    virtual double value() { return val_; }
    virtual ScalarPtr derivative(ScalarPtr other) { return ScalarPtr(new JustScalar(other.get() == this ? 1.0 : 0.0)); }
    virtual VectorPtr derivative(VectorPtr other) { return VectorPtr(new ZeroVector()); }

private:
    double val_;
};

class ScalarAdd : public ScalarScalar2Scalar
{
public:
    ScalarAdd(ScalarPtr a = ScalarPtr(), ScalarPtr b = ScalarPtr()) : ScalarScalar2Scalar()
    {
        setInputScalar1(a);
        setInputScalar2(b);
    }

    // Scalar functions
    virtual double value()
    {
        check(inN1_, "input number 1 not set");
        check(inN2_, "input number 2 not set");

        return inN1_->value() + inN2_->value();
    }

    virtual ScalarPtr derivative(ScalarPtr other)
    {
        if(other.get() == this)
            return ScalarPtr(new JustScalar(1));

        return ScalarPtr(new ScalarAdd(inN1_->derivative(other), inN2_->derivative(other)));
    }

    virtual VectorPtr derivative(VectorPtr other)
    {
        //TBD
        return VectorPtr();
    }
};

inline
ScalarPtr operator+(ScalarPtr a, ScalarPtr b) { return ScalarPtr(new ScalarAdd(a, b)); }

class ScalarSubtract : public ScalarScalar2Scalar
{
public:
    ScalarSubtract(ScalarPtr a = ScalarPtr(), ScalarPtr b = ScalarPtr()) : ScalarScalar2Scalar()
    {
        setInputScalar1(a);
        setInputScalar2(b);
    }

    // Scalar functions
    virtual double value()
    {
        check(inN1_, "input number 1 not set");
        check(inN2_, "input number 2 not set");

        return inN1_->value() - inN2_->value();
    }

    virtual ScalarPtr derivative(ScalarPtr other)
    {
        if(other.get() == this)
            return ScalarPtr(new JustScalar(1));

        return ScalarPtr(new ScalarSubtract(inN1_->derivative(other), inN2_->derivative(other)));
    }

    virtual VectorPtr derivative(VectorPtr other)
    {
        //TBD
        return VectorPtr();
    }
};

inline
ScalarPtr operator-(ScalarPtr a, ScalarPtr b) { return ScalarPtr(new ScalarSubtract(a, b)); }

class ScalarMultiply : public ScalarScalar2Scalar
{
public:
    ScalarMultiply(ScalarPtr a = ScalarPtr(), ScalarPtr b = ScalarPtr()) : ScalarScalar2Scalar()
    {
        setInputScalar1(a);
        setInputScalar2(b);
    }

    // Scalar functions
    virtual double value()
    {
        check(inN1_, "input number 1 not set");
        check(inN2_, "input number 2 not set");

        return inN1_->value() * inN2_->value();
    }

    virtual ScalarPtr derivative(ScalarPtr other)
    {
        if(other.get() == this)
            return ScalarPtr(new JustScalar(1));

        return ScalarPtr(new ScalarMultiply(inN1_->derivative(other), inN2_)) + ScalarPtr(new ScalarMultiply(inN1_, inN2_->derivative(other)));
    }

    virtual VectorPtr derivative(VectorPtr other)
    {
        //TBD
        return VectorPtr();
    }
};

inline
ScalarPtr operator*(ScalarPtr a, ScalarPtr b) { return ScalarPtr(new ScalarMultiply(a, b)); }

class ConstTimesScalar : public Scalar2Scalar
{
public:
    ConstTimesScalar(double c = 1.0, ScalarPtr x = ScalarPtr()) : Scalar2Scalar(), c_(c)
    {
        setInputScalar(x);
    }

    // Scalar functions
    virtual double value()
    {
        check(inN_, "input number not set");
        return c_ * inN_->value();
    }

    virtual ScalarPtr derivative(ScalarPtr other)
    {
        if(other.get() == this)
            return ScalarPtr(new JustScalar(1));

        return ScalarPtr(new ConstTimesScalar(c_, inN_->derivative(other)));
    }

    virtual VectorPtr derivative(VectorPtr other)
    {
        //TBD
        return VectorPtr();
    }

private:
    const double c_;
};

ScalarPtr operator*(double a, ScalarPtr b) { return ScalarPtr(new ConstTimesScalar(a, b)); }

class ScalarPow : public Scalar2Scalar
{
public:
    ScalarPow(ScalarPtr x = ScalarPtr(), double p = 1) : Scalar2Scalar(), p_(p)
    {
        setInputScalar(x);
    }

    // Scalar functions
    virtual double value()
    {
        check(inN_, "input not set");
        return std::pow(inN_->value(), p_);
    }
    
    virtual ScalarPtr derivative(ScalarPtr other)
    {
        if(other.get() == this)
            return ScalarPtr(new JustScalar(1));

        return p_ * ScalarPtr(new ScalarPow(inN_, p_ - 1)) * inN_->derivative(other);
    }

    virtual VectorPtr derivative(VectorPtr other)
    {
        //TBD
        return VectorPtr();
    }

private:
    const double p_;
};

inline
ScalarPtr pow(ScalarPtr x, double p) { return ScalarPtr(new ScalarPow(x, p)); }

class ScalarExp : public Scalar2Scalar
{
public:
    ScalarExp(ScalarPtr x = ScalarPtr()) : Scalar2Scalar()
    {
        setInputScalar(x);
    }

    // Scalar functions
    virtual double value()
    {
        check(inN_, "input not set");
        return std::exp(inN_->value());
    }
    
    virtual ScalarPtr derivative(ScalarPtr other)
    {
        if(other.get() == this)
            return ScalarPtr(new JustScalar(1));

        return ScalarPtr(new ScalarExp(inN_)) * inN_->derivative(other);
    }

    virtual VectorPtr derivative(VectorPtr other)
    {
        //TBD
        return VectorPtr();
    }
};

inline
ScalarPtr exp(ScalarPtr x) { return ScalarPtr(new ScalarExp(x)); }

} // namespace Expressions

#endif


#ifndef COSMO_PP_EXPRESSIONS_HPP
#define COSMO_PP_EXPRESSIONS_HPP

#include <memory>
#include <vector>
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

class VectorDevice
{
public:
    virtual ~VectorDevice() {}

    virtual int create() = 0;
    virtual void remove(int v) = 0;

    virtual void copy(int from, int to) = 0;
    virtual void setToZero(int v) = 0;
    virtual void multiplyBy(int v, double c) = 0;
    virtual void add(int a, int b) = 0; // add b to a
    virtual void swap(int a, int b) = 0;
    virtual void elementwiseMultiply(int a, int b) = 0; // a * b stored in a
    virtual double sum(int v) = 0;
    virtual double dotProduct(int a, int b) = 0;
};

class Vector
{
public:
    virtual ~Vector() {}

    virtual VectorDevice& device() = 0;
    virtual int id() = 0;

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

    void setInputVector(VectorPtr v) { inV_ = v; }

protected:
    VectorPtr inV_;
};

class InVectorVector
{
public:
    virtual ~InVectorVector() {}

    void setInputVector1(VectorPtr v) { inV1_ = v; }
    void setInputVector2(VectorPtr v) { inV2_ = v; }

protected:
    VectorPtr inV1_;
    VectorPtr inV2_;
};

class InScalarVector : public InScalar, public InVector
{
public:
    virtual ~InScalarVector() {}
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

    // Scalar functions
    virtual double value() = 0;
    virtual ScalarPtr derivative(ScalarPtr other) = 0;
    virtual VectorPtr derivative(VectorPtr other) = 0;
};

class Vector2Vector : public Vector, public InVector
{
public:
    virtual ~Vector2Vector() {}

    // Vector functions
    virtual VectorDevice& device() = 0;
    virtual int id() = 0;
    virtual VectorPtr derivative(ScalarPtr other) = 0;
    virtual Vector2VectorPtr derivative(VectorPtr other) = 0;
};

class VectorVector2Vector : public Vector, public InVectorVector
{
public:
    virtual ~VectorVector2Vector() {}

    // Vector functions
    virtual VectorDevice& device() = 0;
    virtual int id() = 0;
    virtual VectorPtr derivative(ScalarPtr other) = 0;
    virtual Vector2VectorPtr derivative(VectorPtr other) = 0;
};

class ScalarVector2Vector : public Vector, public InScalarVector
{
public:
    virtual ~ScalarVector2Vector() {}

    // Vector functions
    virtual VectorDevice& device() = 0;
    virtual int id() = 0;
    virtual VectorPtr derivative(ScalarPtr other) = 0;
    virtual Vector2VectorPtr derivative(VectorPtr other) = 0;
};

class JustScalar : public Scalar
{
public:
    JustScalar(double val = 0) : val_(val) {}
    virtual ~JustScalar() {}

    void set(double val) { val_ = val; }
    virtual double value() { return val_; }
    virtual ScalarPtr derivative(ScalarPtr other) { return ScalarPtr(new JustScalar(other.get() == this ? 1.0 : 0.0)); }
    virtual VectorPtr derivative(VectorPtr other);

private:
    double val_;
};

class JustVector : Vector
{
public:
    JustVector(VectorDevice& d) : d_(d), v_(d.create()) {}
    virtual ~JustVector() { d_.remove(v_); }

    virtual VectorDevice& device() { return d_; }
    virtual int id() { return v_; }

    virtual VectorPtr derivative(ScalarPtr other) { return ScalarPtr(new JustScalar(0)); }
    virtual Vector2VectorPtr derivative(VectorPtr other);
    
private:
    VectorDevice& d_;
    const int v_;
};

VectorPtr
JustScalar::derivative(VectorPtr other);
{
    VectorDevice& d = other->device();
    return VectorPtr<new JustVector(d);
}

class IdentityMatrix : public Vector2Vector
{
public:
    IdentityMatrix(VectorDevice& d) : d_(d) {}
    virtual ~IdentityMatrix() {}

    // Vector functions
    virtual VectorDevice& device() { return d_; }
    virtual int id() { return inV_->id(); }
    virtual VectorPtr derivative(ScalarPtr other) { return inV->derivative(other); }
    virtual Vector2VectorPtr derivative(VectorPtr other)
    {
        if(other.get() == this)
            return Vector2VectorPtr(new IdentityMatrix(d_));

        return inV_->derivative(other);
    }
    
private:
    VectorDevice& d_;
};

class ZeroMatrix : public Vector2Vector
{
public:
    ZeroMatrix(VectorDevice& d) : d_(d), v_(d.create()) {}
    virtual ~ZeroMatrix() { d_.remove(v_); }

    // Vector functions
    virtual VectorDevice& device() { return d_; }
    virtual int id() { return v_; }
    virtual VectorPtr derivative(ScalarPtr other) { return ScalarPtr(new JustScalar(0)); }
    virtual Vector2VectorPtr derivative(VectorPtr other)
    {
        if(other.get() == this)
            return Vector2VectorPtr(new IdentityMatrix(d_));

        return Vector2VectorPtr(new ZeroMatrix(d_));
    }
    
private:
    VectorDevice& d_;
    int v_;
};

Vector2VectorPtr JustVector::derivative(VectorPtr other)
{
    if(other.get() == this)
        return Vector2VectorPtr(new IdentityMatrix(d_));
    return Vector2VectorPtr(new ZeroMatrix(d_));
}


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

class BasicVectorDevice : public VectorDevice
{
public:
    BasicVectorDevice(int n) : n_(n) { check(n_ > 0, ""); }
    virtual ~BasicVectorDevice() {}

    virtual int create()
    {
        data_.push_back(std::vector<double>(n_));
        return data_.size() - 1;
    }

    virtual void remove(int v)
    {
        checkElement(v);
        std::vector<double>().swap(data_[v]); // clear with reallocation
    }

    virtual void copy(int from, int to)
    {
        checkElement(from);
        checkElement(to);
        data_[to] = data_[from];
    }

    virtual void setToZero(int v)
    {
        checkElement(v);
        for(double& x : data_[v])
            x = 0;
    }

    virtual void multiplyBy(int v, double c)
    {
        checkElement(v);
        for(double& x : data_[v])
            x *= c;
    }

    virtual void add(int a, int b) // add b to a
    {
        checkElement(a);
        checkElement(b);

        for(int i = 0; i < n_; ++i)
            a[i] += b[i];
    }

    virtual void swap(int a, int b)
    {
        checkElement(a);
        checkElement(b);
        data_[a].swap(data_[b]);
    }

    virtual void elementwiseMultiply(int a, int b) // a * b stored in a
    {
        checkElement(a);
        checkElement(b);

        for(int i = 0; i < n_; ++i)
            a[i] *= b[i];
    }

    virtual double sum(int v)
    {
        checkElement(v);
        double s = 0;
        for(double x : data_[v])
            s += x;

        // do MPI
        
        return s;
    }

    virtual double dotProduct(int a, int b)
    {
        checkElement(a);
        checkElement(b);
        double s = 0;
        for(int i = 0; i < n_; ++i)
            s += data_[a][i] * data_[b][i];

        // do MPI

        return s;
    }

private:
    void checkElement(int i)
    {
        check(i >= 0 && i < data_.size() && data_[i].size() == n_, "");
    }

private:
    const int n_;
    std::vector<std::vector<double>> data_;
};

} // namespace Expressions

#endif


#ifndef COSMO_PP_FUNCTION_HPP
#define COSMO_PP_FUNCTION_HPP

namespace Math
{

/// An abstract template class for a 1 variable function.
template<typename VarType, typename ValType>
class Function
{
public:
    /// The variable type.
	typedef VarType VariableType;

    /// The value type.
	typedef ValType ValueType;
	
public:
    /// Constructor.

    /// Constructor, does nothing.
	Function() {}

    /// Destructor.

    /// Destructor, does nothing.
	virtual ~Function() {}
	
    /// Evaluate the function.

    /// A pure virtual function for evaluating the function.
    /// \param x The argument of the function.
    /// \return The value of the function.
	virtual ValueType evaluate(VariableType x) const = 0;
};
	
/// An abstract template class for a 2 variable function.
template<typename Var1Type, typename Var2Type, typename ValType>
class Function2
{
public:
    /// The variable 1 type.
	typedef Var1Type Variable1Type;

    /// The variable 2 type.
	typedef Var2Type Variable2Type;

    /// The value type.
	typedef ValType ValueType;
	
public:
    /// Constructor.

    /// Constructor, does nothing.
	Function2() {}

    /// Destructor.

    /// Destructor, does nothing.
	virtual ~Function2() {}
	
    /// Evaluate the function.

    /// A pure virtual function for evaluating the function.
    /// \param x1 Argument 1.
    /// \param x2 Argument 2.
    /// \return The value of the function.
	virtual ValueType evaluate(Variable1Type x1, Variable2Type x2) const = 0;
};

/// An abstract template class for a 3 variable function.
template<typename Var1Type, typename Var2Type, typename Var3Type, typename ValType>
class Function3
{
public:
    /// The variable 1 type.
	typedef Var1Type Variable1Type;

    /// The variable 2 type.
	typedef Var2Type Variable2Type;

    /// The variable 3 type.
	typedef Var3Type Variable3Type;

    /// The value type.
	typedef ValType ValueType;
	
public:
    /// Constructor.

    /// Constructor, does nothing.
	Function3() {}

    /// Destructor.

    /// Destructor, does nothing.
	virtual ~Function3() {}
	
    /// Evaluate the function.

    /// A pure virtual function for evaluating the function.
    /// \param x1 Argument 1.
    /// \param x2 Argument 2.
    /// \param x3 Argument 3.
    /// \return The value of the function.
	virtual ValueType evaluate(Variable1Type x1, Variable2Type x2, Variable3Type x3) const = 0;
};
	
	
/// Real 1 variable function.	
typedef Function<double, double> RealFunction;
	
} //namespace Math

#endif

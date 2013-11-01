#ifndef COSMO_CPP_TABLE_FUNCTION_HPP
#define COSMO_CPP_TABLE_FUNCTION_HPP

#include <istream>
#include <map>

#include <function.hpp>
#include <macros.hpp>

namespace Math
{

/// A class for linear interpolation in one dimension. 
/// It's inherited from both std::map and Function. The map functionality can be used to define and access interpolation points. 
/// The class's functionality can be used to calculate linear interpolation between the points.
template<typename VarType, typename ValType, class Compare = std::less<VarType> >
class TableFunction : public Function<VarType, ValType>, public std::map<VarType, ValType, Compare>
{
public:
	typedef Function<VarType, ValType> BaseType;

    /// The variable type.
	typedef typename BaseType::VariableType VariableType;

    /// The value type.
	typedef typename BaseType::ValueType ValueType;
	
	typedef std::map<VarType, ValType, Compare> MapType;

    /// Iterator type.
	typedef typename MapType::iterator iterator;

    /// Constant iterator type.
	typedef typename MapType::const_iterator const_iterator;
	
public:
    /// Constructor.
    /// \param comp The compare (less) predicate, if needed.
	TableFunction(const Compare& comp = Compare()) : MapType(comp) {}

    /// Destructor.
	~TableFunction() {}
	
    /// Evaluate the linear interpolation.
    /// \param x The argument. Must be between the lowest and highest points defined.
    /// \return The value of the interpolation.
	ValueType evaluate(VariableType x) const;
};
	
template<typename VarType, typename ValType, class Compare>
typename TableFunction<VarType, ValType, Compare>::ValueType TableFunction<VarType, ValType, Compare>::evaluate(VariableType x) const
{
	check(!MapType::empty(), "The map is empty!");
	const_iterator it = MapType::lower_bound(x);
	check(it != MapType::end(), "Element is outside the range!");
	const VariableType x2 = (*it).first;
	const ValueType y2 = (*it).second;
	check(!MapType::key_comp()(x2, x), "");
	if(!MapType::key_comp()(x, x2))
		return y2;
	check(it != MapType::begin(), "Element is outside the range!");
	--it;
	const VariableType x1 = (*it).first;
	const ValueType y1 = (*it).second;
	check(MapType::key_comp()(x1, x), "");
	
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}
	
/// A class for linear interpolation in two dimension. 
/// It's inherited from both std::map and Function. The map functionality can be used to define and access interpolation points. 
/// The class's functionality can be used to calculate linear interpolation between the points.
/// The grid defined needs to be rectangular, i.e. if the points (x1, y1) and (x2, y2) are defined then so should be (x1, y2) and (x2, y1).
template<typename Var1Type, typename Var2Type, typename ValType, class Compare1 = std::less<Var1Type>, class Compare2 = std::less<Var2Type> >
class TableFunction2 : public Function2<Var1Type, Var2Type, ValType>, public std::map<Var1Type, TableFunction<Var2Type, ValType, Compare2>, Compare1>
{
public:
	typedef Function2<Var1Type, Var2Type, ValType> BaseType;

    /// The first variable type.
	typedef typename BaseType::Variable1Type Variable1Type;

    /// The second variable type.
	typedef typename BaseType::Variable2Type Variable2Type;

    /// The value type.
	typedef typename BaseType::ValueType ValueType;
	
	typedef TableFunction<Var2Type, ValType, Compare2> TF1Type;
	typedef std::map<Var1Type, TF1Type, Compare1> MapType;
    
    /// Iterator type.
	typedef typename MapType::iterator iterator;

    /// Constant iterator type.
	typedef typename MapType::const_iterator const_iterator;
	
public:
    /// Constructor.
    /// \param comp1 The predicate (less) for variable 1, if needed.
    /// \param comp2 The predicate (less) for variable 2, if needed.
	TableFunction2(const Compare1& comp1 = Compare1(), const Compare2& comp2 = Compare2()) : MapType(comp1), default_(comp2) {}

    /// Destructor.
	~TableFunction2() {}
	
    /// Access to the one dimensional table function for the second variable, given the fixed value of the first variable.
    /// \param var The first variable.
    /// \return A reference to the table function for the second variable.
	TF1Type& operator[](const Variable1Type& var)
	{
		iterator it = MapType::find(var);
		if(it != MapType::end())
			return it->second;
		
		std::pair<iterator, bool> r = MapType::insert(std::make_pair(var, default_));
		check(r.second, "");
		return (r.first)->second;
	}
	
    /// Evaluate the linear interpolation.
    /// \param x1 The first variable.
    /// \param x2 The second variable.
    /// \return The value of the linear interpolation.
	ValueType evaluate(Variable1Type x1, Variable2Type x2) const;
	
private:
	TF1Type default_;
};

template<typename Var1Type, typename Var2Type, typename ValType, class Compare1, class Compare2>
typename TableFunction2<Var1Type, Var2Type, ValType, Compare1, Compare2>::ValueType
	TableFunction2<Var1Type, Var2Type, ValType, Compare1, Compare2>::evaluate(Variable1Type x1, Variable2Type x2) const
{
	check(!MapType::empty(), "The map is empty!");
	const_iterator it = lower_bound(x1);
	check(it != MapType::end(), "Element is outside the range!");
	const Variable1Type x12 = (*it).first;
	const ValueType y2 = (*it).second.evaluate(x2);
	check(!MapType::key_comp()(x12, x1), "");
	if(!MapType::key_comp()(x1, x12))
		return y2;
	check(it != MapType::begin(), "Element is outside the range!");
	--it;
	const Variable1Type x11 = (*it).first;
	const ValueType y1 = (*it).second.evaluate(x2);
	check(MapType::key_comp()(x11, x1), "");
	
	return y1 + (x1 - x11) * (y2 - y1) / (x12 - x11);
}



/// A class for linear interpolation in three dimension. 
/// It's inherited from both std::map and Function. The map functionality can be used to define and access interpolation points. 
/// The class's functionality can be used to calculate linear interpolation between the points.
/// The grid defined needs to be rectangular, i.e. if the points (x1, y1, z1) and (x2, y2, z2) are defined then so should be 
/// (x1, y1, z2), (x1, y2, z1), (x1, y2, z2), (x2, y1, z1), (x2, y1, z2), and (x2, y2, z1).
template<typename Var1Type, typename Var2Type, typename Var3Type, typename ValType, 
	class Compare1 = std::less<Var1Type>, class Compare2 = std::less<Var2Type>, class Compare3 = std::less<Var3Type> >
class TableFunction3 : public Function3<Var1Type, Var2Type, Var3Type, ValType>, 
	public std::map<Var1Type, TableFunction2<Var2Type, Var3Type, ValType, Compare2, Compare3>, Compare1>
{
public:
	typedef Function3<Var1Type, Var2Type, Var3Type, ValType> BaseType;

    /// The first variable type.
	typedef typename BaseType::Variable1Type Variable1Type;

    /// The second variable type.
	typedef typename BaseType::Variable2Type Variable2Type;

    /// The third variable type.
	typedef typename BaseType::Variable3Type Variable3Type;

    /// The value type.
	typedef typename BaseType::ValueType ValueType;
	
	typedef TableFunction2<Var2Type, Var3Type, ValType, Compare2, Compare3> TF2Type;
	typedef std::map<Var1Type, TF2Type, Compare1> MapType;

    /// Iterator type.
	typedef typename MapType::iterator iterator;

    /// Constant iterator type.
	typedef typename MapType::const_iterator const_iterator;
	
public:
    /// Constructor.
    /// \param comp1 The predicate (less) for the first variable, if needed.
    /// \param comp2 The predicate (less) for the second variable, if needed.
    /// \param comp3 The predicate (less) for the third variable, if needed.
	TableFunction3(const Compare1& comp1 = Compare1(), const Compare2& comp2 = Compare2(), const Compare3& comp3 = Compare3())
	: MapType(comp1), default_(comp2, comp3) {}

    /// Destructor.
	~TableFunction3() {}
	
    /// Access to the two dimensional table function for the second and third variables, given the fixed value of the first variable.
    /// \param var The first variable.
    /// \return A reference to the table function for the second and third variables.
	TF2Type& operator[](const Variable1Type& var)
	{
		iterator it = MapType::find(var);
		if(it != MapType::end())
			return it->second;
		
		std::pair<iterator, bool> r = MapType::insert(std::make_pair(var, default_));
		check(r.second, "");
		return (r.first)->second;
	}
	
    /// Evaluate the linear interpolation.
    /// \param x1 The first variable.
    /// \param x2 The second variable.
    /// \param x3 The third variable.
    /// \return The value of the linear interpolation.
	ValueType evaluate(Variable1Type x1, Variable2Type x2, Variable3Type x3) const;
	
private:
	TF2Type default_;
};

template<typename Var1Type, typename Var2Type, typename Var3Type, typename ValType, class Compare1, class Compare2, class Compare3>
typename TableFunction3<Var1Type, Var2Type, Var3Type, ValType, Compare1, Compare2, Compare3>::ValueType
	TableFunction3<Var1Type, Var2Type, Var3Type, ValType, Compare1, Compare2, Compare3>::evaluate(Variable1Type x1, Variable2Type x2, Variable3Type x3) const
{
	check(!MapType::empty(), "The map is empty!");
	const_iterator it = lower_bound(x1);
	check(it != MapType::end(), "Element is outside the range!");
	const Variable1Type x12 = (*it).first;
	const ValueType y2 = (*it).second.evaluate(x2, x3);
	check(!MapType::key_comp()(x12, x1), "");
	if(!MapType::key_comp()(x1, x12))
		return y2;
	check(it != MapType::begin(), "Element is outside the range!");
	--it;
	const Variable1Type x11 = (*it).first;
	const ValueType y1 = (*it).second.evaluate(x2, x3);
	check(MapType::key_comp()(x11, x1), "");
	
	return y1 + (x1 - x11) * (y2 - y1) / (x12 - x11);
}
	
} //namespace Math

#endif

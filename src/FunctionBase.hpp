/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef FUNCTIONBASE_HPP
#define FUNCTIONBASE_HPP

namespace Galerkin
{

/*!
 * @brief Functions contains functionality for **scalar** functions
 * 
 * The utilities in this namespace exist to do calculus and arithmetic with
 * functions from R^n -> R. It contains a base class that serves to provide
 * fallback methods for summation and products of functions, and to compute their
 * partial derivatives.
 */
namespace Functions
{

template <class Derived>
class FunctionBase;

template <class F1, class F2>
class FunctionSum : public FunctionBase<FunctionSum<F1, F2>>
{
public:
    constexpr FunctionSum(const F1 &a, const F2 &b) : m_f1(a), m_f2(b) {}

    template <auto I>
    constexpr auto partial() const noexcept
    {
        return partial<I>(m_f1) + partial<I>(m_f2);
    }

    template <class... T>
    constexpr auto operator()(const T &... x) const noexcept
    {
        return m_f1(x...) + m_f2(x...);
    }

private:
    F1 m_f1;
    F2 m_f2;
};

template <class F1, class F2>
class FunctionProduct : public FunctionBase<FunctionProduct<F1, F2>>
{
public:
    constexpr FunctionProduct(const F1 &a, const F2 &b) : m_f1(a), m_f2(b) {}

    template <auto I>
    constexpr auto partial() const noexcept
    {
        return partial<I>(m_f1) * m_f2 + m_f1 * partial<I>(m_f2);
    }

    template <class... T>
    constexpr auto operator()(const T &... x) const noexcept
    {
        return m_f1(x...) * m_f2(x...);
    }

private:
    F1 m_f1;
    F2 m_f2;
};

template <class Derived>
class FunctionBase
{
public:
    template <class OtherDerived>
    constexpr auto operator+(const FunctionBase<OtherDerived> &other) const noexcept
    {
        return FunctionSum(*this, other);
    }

    template <class OtherDerived>
    constexpr auto operator*(const FunctionBase<OtherDerived> &other) const noexcept
    {
        return FunctionProduct(*this, other);
    }

private:
    FunctionBase() = default;
    friend Derived;
};

} // namespace Functions

} // namespace Galerkin

#endif /* FUNCTIONBASE_HPP */
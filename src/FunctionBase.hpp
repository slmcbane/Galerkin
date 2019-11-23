/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef FUNCTIONBASE_HPP
#define FUNCTIONBASE_HPP

#include "utils.hpp"
#include "Rationals.hpp"

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
        return Galerkin::partial<I>(m_f1) + Galerkin::partial<I>(m_f2);
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
        return Galerkin::partial<I>(m_f1) * m_f2 + m_f1 * Galerkin::partial<I>(m_f2);
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

template <class F>
class FunctionNegation : public FunctionBase<FunctionNegation<F>>
{
public:
    constexpr FunctionNegation(const F& f) : m_f(f) {}

    template <auto I>
    constexpr auto partial() const noexcept
    {
        return FunctionNegation(Galerkin::partial<I>(m_f));
    }

    template <class... T>
    constexpr auto operator()(const T&... x) const noexcept
    {
        return -m_f(x...);
    }

private:
    F m_f;
};

template <class F1, class F2>
class FunctionQuotient : public FunctionBase<FunctionQuotient<F1, F2>>
{
public:
    constexpr FunctionQuotient(const F1 &a, const F2 &b) : m_f1(a), m_f2(b) {}

    template <auto I>
    constexpr auto partial() const noexcept
    {
        return (m_f2 * Galerkin::partial<I>(m_f1) -
            m_f1 * Galerkin::partial<I>(m_f2)) / (m_f2 * m_f2);
    }

    template <class... T>
    constexpr auto operator()(const T &... x) const noexcept
    {
        return m_f1(x...) / m_f2(x...);
    }

private:
    F1 m_f1;
    F2 m_f2;
};

template <class Derived>
class FunctionBase
{
public:
    template <class Other>
    constexpr auto operator+(const Other &other) const noexcept
    {
        return FunctionSum(static_cast<const Derived&>(*this), other);
    }

    template <class Other>
    constexpr auto operator-(const Other &other) const noexcept
    {
        return *this + FunctionNegation(other);
    }

    constexpr auto operator-() const noexcept
    {
        return FunctionNegation(static_cast<const Derived&>(*this));
    }

    template <class Other>
    constexpr auto operator*(const Other &other) const noexcept
    {
        return FunctionProduct(static_cast<const Derived&>(*this), other);
    }

    template <class Other>
    constexpr auto operator/(const Other &other) const noexcept
    {
        return FunctionQuotient(static_cast<const Derived&>(*this), other);
    }

    template <class... T>
    constexpr auto operator()(const T &... x) const noexcept
    {
        return static_cast<const Derived*>(*this)->operator()(x...);
    }

private:
    FunctionBase() = default;
    friend Derived;
};

template <class T>
class ConstantFunction : public FunctionBase<ConstantFunction<T>>
{
public:
    constexpr ConstantFunction(T v) noexcept : m_val(v) {}

    template <class... Args>
    constexpr auto operator()([[maybe_unused]] const Args&... args) const noexcept
    {
        return m_val;
    }

    template <auto I>
    constexpr auto partial() const noexcept
    {
        return ConstantFunction(Rationals::rational<0>);
    }
private:
    T m_val;
};

/********************************************************************************
 * Tests for ConstantFunction
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Functions] Test ConstantFunction")
{
    constexpr auto fn = ConstantFunction(2.0);
    REQUIRE(fn() == 2.0);
    REQUIRE(fn(1, 2, 3) == 2.0);
    REQUIRE(fn(std::tuple('a', 'b', 'c')) == 2.0);

    REQUIRE(partial<0>(fn)(1) == Rationals::rational<0>);
    REQUIRE(partial<10>(fn)(std::tuple(0, 0)) == Rationals::rational<0>);
}

#endif // DOCTEST_LIBRARY_INCLUDED
/********************************************************************************
 * End test block.
 *******************************************************************************/

} // namespace Functions

} // namespace Galerkin

#endif /* FUNCTIONBASE_HPP */
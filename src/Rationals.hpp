/*
 * Copyright (c) 2021, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef RATIONALS_HPP
#define RATIONALS_HPP

/*!
 * @file Rationals.hpp
 * @brief Compile-time computations with rational numbers.
 * @author Sean McBane <sean.mcbane@protonmail.com>
 */

#include <cstdint>
#include <type_traits>

#include "utils.hpp"

namespace Galerkin
{

/// The `Rationals` namespace encapsulates all of the functionality for rational numbers
namespace Rationals
{

/*!
 * @brief The type used for the numerator in a rational. On gcc or clang this
 * could be made a int128_t, but in practice compiler limits will be reached
 * before overflow, anyway. Should be a signed type.
 */
typedef int64_t rational_num_type;

/// The type used for the denominator in a rational. Should be an unsigned type.
typedef uint64_t rational_den_type;

/*!
 * @brief Utility; find greatest common denominator of `a` and `b`. This will hit
 * `constexpr` evaluation limits when `a` or `b` becomes large relative to the
 * other.
 */
constexpr auto gcd(rational_num_type a, rational_den_type b)
{
    if (a < 0)
    {
        a = -a;
    }

    auto x = static_cast<rational_den_type>(a);
    auto y = x > b ? b : x;
    x = x > b ? x : b;

    if (y == 0)
    {
        return x;
    }

    while (x != y)
    {
        x = x - y;
        if (y > x)
        {
            auto tmp = y;
            y = x;
            x = tmp;
        }
    }

    return x;
}

namespace detail
{

/// Reduce a rational number to its simplest representation.
constexpr auto reduce_rational(rational_num_type n, rational_den_type d)
{
    auto div = gcd(n, d);
    return std::tuple(n / static_cast<rational_num_type>(div), d / div);
}

} // namespace detail

/*!
 * @brief Type representing a rational number at compile time.
 *
 * `Rational` is used to do exact calculations of derivatives, Jacobians, etc.
 * by manipulating rational numbers at compile time. The normal arithmetic
 * operators are all defined for it, and it decays to a `double` as appropriate.
 * Rather than use this class template, however, you should use the `rational`
 * template constant from this header; instantiating `rational<N, D>`
 * automatically reduces the resulting fraction to its simplest form.
 */
class Rational
{
  public:
    /// Get the numerator
    constexpr auto num() const noexcept { return m_num; }
    /// Get the denominator
    constexpr auto den() const noexcept { return m_den; }

    /*! @brief Convert the number to a double precision float, for when you have to
     * do inexact floating point operations or run-time computation :(.
     */
    constexpr operator double() const noexcept { return static_cast<double>(m_num) / m_den; }

    template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
    constexpr Rational(T x) : m_num(x), m_den{1}
    {
    }

    template <
        class T1, class T2, class = std::enable_if_t<std::is_integral_v<T1> && std::is_integral_v<T2>>>
    constexpr Rational(T1 x, T2 y) : m_num(x), m_den{1}
    {
        if (y < 0)
        {
            m_num *= -1;
            m_den = -y;
        }
        else
        {
            m_den = y;
        }
        auto reduced = detail::reduce_rational(m_num, m_den);
        m_num = std::get<0>(reduced);
        m_den = std::get<1>(reduced);
    }

  private:
    rational_num_type m_num;
    rational_den_type m_den;
};

/*!
 * @brief Template constant for a compile-time rational.
 *
 * This constant evaluates to a `Rational` reduced to its lowest terms. Default
 * denominator is 1, so that an integer can be constructed by `rational<n>`.
 *
 * @tparam N The numerator
 * @tparam D The denominator. Default value: 1
 */
template <rational_num_type N, rational_den_type D = 1>
constexpr auto rational = Rational(N, D);

constexpr bool operator==(const Rational &x, const Rational &y)
{
    return (x.num() == y.num() && x.den() == y.den());
}

/********************************************************************************
 * Tests of rational construction
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Rationals] Testing construction of rationals")
{
    REQUIRE(rational<1> == rational<2, 2>);
    REQUIRE(rational<-2, 2> == rational<-42, 42>);
    REQUIRE(rational<2, 4> == rational<1, 2>);
    REQUIRE(rational<4, 2> == rational<2, 1>);
    REQUIRE(rational<-1, 3> == rational<-6, 18>);
    // This should trigger a static assert
    // REQUIRE(rational<1, 0> == rational<0, 1>);
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

constexpr auto operator+(const Rational &x, const Rational &y)
{
    auto lcm = x.den() * y.den() / gcd(x.den(), y.den());
    auto mult1 = lcm / x.den();
    auto mult2 = lcm / y.den();

    return Rational(
        x.num() * static_cast<rational_num_type>(mult1) + y.num() * static_cast<rational_num_type>(mult2),
        x.den() * mult1);
}

/********************************************************************************
 * Tests of rational addition
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Rationals] Testing rational addition")
{
    REQUIRE(rational<1> + rational<2> == rational<3>);
    REQUIRE(rational<1, 2> + rational<1, 3> == rational<5, 6>);
    REQUIRE(rational<5, 6> + rational<1, 6> == rational<1>);
    REQUIRE(rational<5, 8> + rational<22, 16> == rational<2>);
    REQUIRE(rational<1, 3> + rational<3, 1> == rational<10, 3>);
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

constexpr auto operator-(const Rational &x) noexcept
{
    return Rational(-x.num(), x.den());
}

constexpr auto operator-(const Rational &x, const Rational &y)
{
    return x + -y;
}

/********************************************************************************
 * Tests of rational subtraction and negation.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Rationals] Testing rational subtraction")
{
    REQUIRE(rational<1> - rational<2> == rational<-1>);
    REQUIRE(rational<1, 2> - rational<1, 3> == rational<1, 6>);
    REQUIRE(rational<5, 6> - rational<1, 6> == rational<2, 3>);
    REQUIRE(rational<5, 8> - rational<22, 16> == rational<-3, 4>);
    REQUIRE(rational<1, 3> - rational<3, 1> == rational<-8, 3>);
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

constexpr auto operator*(const Rational &x, const Rational &y)
{
    return Rational(x.num() * y.num(), x.den() * y.den());
}

constexpr auto operator/(const Rational &x, const Rational &y)
{
    auto num = x.num() * static_cast<rational_num_type>(y.den());
    if (y.num() < 0)
    {
        return Rational(-num, x.den() * static_cast<rational_den_type>(-y.num()));
    }
    else
    {
        return Rational(num, x.den() * y.num());
    }
}

/// A `Rational` * an `integral_constant` returns a `Rational`.
template <class T>
constexpr std::enable_if_t<std::is_integral_v<T>, Rational>
operator*(const Rational &x, T y) noexcept
{
    return Rational(x.num() * static_cast<rational_num_type>(y), x.den());
}

template <class T>
constexpr std::enable_if_t<std::is_integral_v<T>, Rational>
operator*(T y, const Rational &x) noexcept
{
    return x * y;
}

template <class T, T v>
constexpr Rational operator*(const Rational &x, std::integral_constant<T, v>)
{
    return x * v;
}

template <class T, T v>
constexpr Rational operator*(std::integral_constant<T, v>, const Rational &x)
{
    return x * v;
}

/// A `Rational` / an `integral_constant` returns a `Rational`.
template <class T>
constexpr std::enable_if_t<std::is_integral_v<T>, Rational>
operator/(const Rational &x, T y)
{
    if (y < 0)
    {
        return Rational(-x.num(), x.den() * (-y));
    }
    else
    {
        return Rational(x.num(), x.den() * y);
    }
}

template <class T>
constexpr std::enable_if_t<std::is_integral_v<T>, Rational>
operator/(T y, const Rational &x)
{
    if (x.num() < 0)
    {
        return -Rational(y * static_cast<rational_num_type>(x.den()), -x.num());
    }
    else
    {
        return Rational(y * static_cast<rational_num_type>(x.den()), x.num());
    }
}

template <class T, T v>
constexpr Rational operator/(const Rational &x, std::integral_constant<T, v>)
{
    return x / v;
}

template <class T, T v>
constexpr Rational operator/(std::integral_constant<T, v>, const Rational &x)
{
    return v / x;
}

template <class T>
constexpr std::enable_if_t<std::is_floating_point_v<T>, double>
operator*(const Rational &x, T y)
{
    return static_cast<double>(x) * y;
}

template <class T>
constexpr std::enable_if_t<std::is_floating_point_v<T>, double>
operator/(const Rational &x, T y)
{
    return static_cast<double>(x) / y;
}

/********************************************************************************
 * Tests of rational multiplication and division.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Rationals] Testing rational multiplication and division")
{
    REQUIRE(rational<1, 2> * rational<1, 2> == rational<1, 4>);
    REQUIRE(rational<1, 2> * rational<1, 3> == rational<1, 3> * rational<1, 2>);
    REQUIRE(rational<1, 2> * rational<1, 3> == rational<1, 6>);
    REQUIRE(rational<3, 10> * rational<1, 3> == rational<1, 10>);
    REQUIRE(rational<3, 10> * rational<-1, 3> == rational<-1, 10>);

    REQUIRE(rational<1, 2> / rational<1, 2> == rational<1>);
    REQUIRE(rational<1, 2> / rational<2> == rational<1, 4>);
    REQUIRE(rational<3, 10> / rational<1, 3> == rational<9, 10>);
    REQUIRE(rational<1, 6> / rational<1, 3> == rational<1, 2>);
    REQUIRE(rational<1, 6> / rational<1, 2> == rational<1, 3>);
    REQUIRE(rational<3, 10> / rational<-1, 3> == rational<-9, 10>);

    REQUIRE(rational<3, 10> / std::integral_constant<int, 3>() == rational<1, 10>);
    REQUIRE(rational<1, 6> * std::integral_constant<int, 3>() == rational<1, 2>);
    REQUIRE(rational<1, 6> / std::integral_constant<int, -2>() == rational<-1, 12>);

    REQUIRE(rational<1, 2> * 0.5 == doctest::Approx(0.25));
    REQUIRE(rational<1, 2> / 3 == doctest::Approx(1.0 / 6));
    REQUIRE(std::is_same_v<decltype(rational<1, 2> * std::declval<float>()), double>);
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

} // namespace Rationals

} // namespace Galerkin

#endif /* RATIONALS_HPP */
/*
 * Copyright (c) 2019, Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef RATIONALS_HPP
#define RATIONALS_HPP

#include <cstdint>
#include <type_traits>

#include "utils.hpp"

namespace Galerkin
{

namespace Rationals
{

typedef int64_t rational_num_t;
typedef uint64_t rational_den_t;

template <rational_num_t Num, rational_den_t Den>
struct Rational
{
    static_assert(Den != zero<rational_den_t>());

    static constexpr auto num() { return Num; }
    static constexpr auto den() { return Den; }
};

constexpr auto gcd(rational_num_t a, rational_den_t b)
{
    if (a < 0)
    {
        a = -a;
    }

    auto x = a > b ? a : b;
    auto y = a > b ? b : a;

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

template <rational_num_t N, rational_den_t D>
constexpr auto reduce_rational(Rational<N, D>)
{
    constexpr auto div = gcd(N, D);
    return Rational<N / static_cast<rational_num_t>(div), D / div>();
}

template <rational_num_t N, rational_den_t D = 1>
constexpr auto rational = reduce_rational(Rational<N, D>());

template <auto N1, auto D1, auto N2, auto D2>
constexpr bool operator==(Rational<N1, D1>, Rational<N2, D2>)
{
    return std::is_same_v<decltype(rational<N1, D1>), decltype(rational<N2, D2>)>;
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

template <auto N1, auto D1, auto N2, auto D2>
constexpr auto operator+(Rational<N1, D1>, Rational<N2, D2>)
{
    constexpr auto lcm = D1 * D2 / gcd(D1, D2);
    constexpr auto mult1 = lcm / D1;
    constexpr auto mult2 = lcm / D2;

    return rational<N1*mult1 + N2*mult2, D1*mult1>;
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

template <auto N1, auto D1>
constexpr auto operator-(Rational<N1, D1>)
{
    return rational<-N1, D1>;
}

template <auto N1, auto D1, auto N2, auto D2>
constexpr auto operator-(Rational<N1, D1>, Rational<N2, D2>)
{
    return Rational<N1, D1>() + (-Rational<N2, D2>());
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

template <auto N1, auto D1, auto N2, auto D2>
constexpr auto operator*(Rational<N1, D1>, Rational<N2, D2>)
{
    return rational<N1*N2, D1*D2>;
}

template <auto N1, auto D1, auto N2, auto D2>
constexpr auto operator/(Rational<N1, D1>, Rational<N2, D2>)
{
    constexpr auto num = N1 * static_cast<rational_num_t>(D2);
    if constexpr (N2 < 0)
    {
        return -rational<num, D1 * static_cast<rational_den_t>(-N2)>;
    }
    else
    {
        return rational<num, D1 * static_cast<rational_den_t>(N2)>;
    }
}

template <auto N, auto D, class I, I v>
constexpr auto operator*(Rational<N, D>, std::integral_constant<I, v>)
{
    return rational<N, D> * rational<v>;
}

template <auto N, auto D, class I, I v>
constexpr auto operator/(Rational<N, D>, std::integral_constant<I, v>)
{
    static_assert(v != 0);
    if constexpr (v < 0)
    {
        return -(rational<N, D> / std::integral_constant<I, -v>());
    }
    else
    {
        return rational<N, D> / rational<v>;
    }
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
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/


} // namespace Rationals

} // namespace Galerkin

#endif /* RATIONALS_HPP */
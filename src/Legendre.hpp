#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include "Polynomial.hpp"
#include "Rationals.hpp"

#include <array>
#include <limits>

namespace Galerkin
{

namespace Legendre
{

// This works up to 14th order before hitting a constexpr evaluation limit in
// clang.
template <auto I>
constexpr inline auto polynomial = (Polynomials::make_poly(
                                        std::tuple(Rationals::Rational(2*I-1)),
                                        Polynomials::PowersList<Polynomials::Powers<1>>{}) *
                                        polynomial<I - 1> +
                                    polynomial<I - 2> * Rationals::rational<1 - I>) / Rationals::rational<I>;

template <>
constexpr inline auto polynomial<0> = Polynomials::make_poly(
    std::tuple(Rationals::Rational(1)), Polynomials::PowersList<Polynomials::Powers<0>>{});

template <>
constexpr inline auto polynomial<1> = Polynomials::make_poly(
    std::tuple(Rationals::Rational(1)), Polynomials::PowersList<Polynomials::Powers<1>>{});

/********************************************************************************
 * Test that Legendre polynomials are correct.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using namespace Polynomials;
using namespace Rationals;

TEST_CASE("Test computed Legendre polynomials")
{
    SUBCASE("Order 0")
    {
        constexpr auto poly = polynomial<0>;
        REQUIRE(poly.coeffs()[0] == rational<1>);
        static_assert(std::is_same_v<std::decay_t<decltype(poly)>, Polynomial<Rational, Powers<0>>>);
    }

    SUBCASE("Order 1")
    {
        constexpr auto poly = polynomial<1>;
        REQUIRE(poly.coeffs()[0] == rational<1>);
        static_assert(std::is_same_v<std::decay_t<decltype(poly)>, Polynomial<Rational, Powers<1>>>);
    }

    SUBCASE("Order 2")
    {
        constexpr auto poly = polynomial<2>;
        REQUIRE(poly.coeffs()[0] == -rational<1, 2>);
        REQUIRE(poly.coeffs()[1] == rational<3, 2>);
        static_assert(
            std::is_same_v<std::decay_t<decltype(poly)>, Polynomial<Rationals::Rational, Powers<0>, Powers<2>>>);
    }

    SUBCASE("Order 10")
    {
        constexpr auto poly = polynomial<10>;
        REQUIRE(poly.coeffs()[0] == -rational<63, 256>);
        REQUIRE(poly.coeffs()[1] == rational<3465, 256>);
        REQUIRE(poly.coeffs()[2] == -rational<30030, 256>);
        REQUIRE(poly.coeffs()[3] == rational<90090, 256>);
        REQUIRE(poly.coeffs()[4] == -rational<109395, 256>);
        REQUIRE(poly.coeffs()[5] == rational<46189, 256>);
        static_assert(
            std::is_same_v<
                std::decay_t<decltype(poly)>,
                Polynomial<Rational, Powers<0>, Powers<2>, Powers<4>, Powers<6>, Powers<8>, Powers<10>>>);
    }
} // TEST_CASE

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

template <class T, class P>
constexpr auto interval_root(P poly, T low, T high)
{
    auto pprime = poly.template partial<0>();
    T mid = (high + low) / 2;
    constexpr auto abs = [](auto x) { return x < 0 ? -x : x; };
    while (poly(mid) != zero<T>)
    {
        auto delta = -poly(mid) / pprime(mid);
        if (abs(delta / mid) < 4 * std::numeric_limits<T>::epsilon())
        {
            break;
        }
        mid += delta;
    }
    return mid;
}

template <class T, auto Order>
constexpr auto all_roots()
{
    if constexpr (Order == 1)
    {
        return std::array<T, 1>{0};
    }
    else
    {
        auto extrema = all_roots<T, Order - 1>();
        std::array<T, Order> my_roots{0};
        my_roots[0] = interval_root(polynomial<Order>, -one<T>, extrema[0]);
        for (int i = 1; i < Order - 1; ++i)
        {
            my_roots[i] = interval_root(polynomial<Order>, extrema[i - 1], extrema[i]);
        }
        my_roots[Order - 1] = interval_root(polynomial<Order>, extrema[Order - 2], one<T>);
        return my_roots;
    }
}

// This works up to order 6 before hitting a constexpr evaluation limit for my
// tested version of clang.
template <class T, auto Order>
constexpr auto roots = all_roots<T, Order>();

/********************************************************************************
 * Test rootfinding of Legendre polynomials.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

#include <cmath>

TEST_CASE("Find roots of Legendre polynomials")
{
    REQUIRE(roots<double, 1>[0] == doctest::Approx(0.0));

    REQUIRE(roots<double, 2>[0] == doctest::Approx(-1 / std::sqrt(3)));
    REQUIRE(roots<double, 2>[1] == doctest::Approx(1 / std::sqrt(3)));

    REQUIRE(roots<double, 3>[0] == doctest::Approx(-std::sqrt(3.0 / 5)));
    REQUIRE(roots<double, 3>[1] == doctest::Approx(0.0));
    REQUIRE(roots<double, 3>[2] == doctest::Approx(std::sqrt(3.0 / 5)));

    REQUIRE(roots<double, 5>[0] == doctest::Approx(-std::sqrt(5 + 2 * std::sqrt(10.0 / 7)) / 3));
    REQUIRE(roots<double, 5>[1] == doctest::Approx(-std::sqrt(5 - 2 * std::sqrt(10.0 / 7)) / 3));
    REQUIRE(roots<double, 5>[2] == doctest::Approx(0.0));
    REQUIRE(roots<double, 5>[3] == doctest::Approx(std::sqrt(5 - 2 * std::sqrt(10.0 / 7)) / 3));
    REQUIRE(roots<double, 5>[4] == doctest::Approx(std::sqrt(5 + 2 * std::sqrt(10.0 / 7)) / 3));
}

#endif /* DOCTEST_LIBRARY_INCLUDED *

/********************************************************************************
 *******************************************************************************/
} /* namespace Legendre */

} /* namespace Galerkin */

#endif /* LEGENDRE_HPP */
/*
 * Copyright (c) 2019, Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <array>
#include "utils.hpp"
#include "Legendre.hpp"

namespace Galerkin
{

namespace Quadrature
{

template <class T, auto N>
struct Rule
{
    const std::array<T, N> points;
    const std::array<T, N> weights;
    const std::array<T, 2> interval;
};

template <class F, class T, auto N>
constexpr auto integrate(const F &f, Rule<T, N> rule) noexcept
{
    T x = zero<T>;
    for (int i = 0; i < N; ++i)
    {
        x += f(rule.points[i]) * rule.weights[i];
    }
    return x;
}

template <class T, auto N>
constexpr auto legendre_weights() noexcept
{
    std::array<T, N> weights{};
    constexpr std::array<T, N> roots = Legendre::roots<T, N>;
    constexpr auto pprime = partial<0>(Legendre::polynomial<N>);
    int i = 0;
    for (T x : roots)
    {
        weights[i++] = 2 / ((1 - x * x) * pprime(x) * pprime(x));
    }
    return weights;
}

/********************************************************************************
 * Test that computed quadrature points are correct.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

#include <cmath>

TEST_CASE("Test compute points for Gauss-Legendre quadrature")
{
    REQUIRE(legendre_weights<double, 1>()[0] == doctest::Approx(2.0));
    REQUIRE(legendre_weights<double, 2>()[0] == doctest::Approx(1.0));
    REQUIRE(legendre_weights<double, 2>()[1] == doctest::Approx(1.0));
    REQUIRE(legendre_weights<double, 3>()[0] == doctest::Approx(5.0 / 9));
    REQUIRE(legendre_weights<double, 3>()[1] == doctest::Approx(8.0 / 9));
    REQUIRE(legendre_weights<double, 3>()[2] == doctest::Approx(5.0 / 9));
    REQUIRE(legendre_weights<double, 4>()[0] == doctest::Approx((18 - std::sqrt(30)) / 36));
    REQUIRE(legendre_weights<double, 4>()[1] == doctest::Approx((18 + std::sqrt(30)) / 36));
    REQUIRE(legendre_weights<double, 4>()[2] == doctest::Approx((18 + std::sqrt(30)) / 36));
    REQUIRE(legendre_weights<double, 4>()[3] == doctest::Approx((18 - std::sqrt(30)) / 36));
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 * End test of computed quadrature points.
 *******************************************************************************/

template <class T, auto N>
constexpr Rule<T, N> legendre_rule = Rule<T, N> { Legendre::roots<T, N>,
                                                  legendre_weights<T, N>(),
                                                  std::array<T, 2>{-one<T>, one<T>} };

/********************************************************************************
 * Test that an integral is computed accurately.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Test that integrals are computed exactly")
{
    SUBCASE("Check the integral of an order 1 polynomial with 1 quadrature point")
    {
        REQUIRE(
            integrate( [](double x) { return 3*x + 2; }, legendre_rule<double, 1> )
                == doctest::Approx(3.0 / 2 + 2 - (3.0 / 2 - 2))
        );
    }

    SUBCASE("Check the integral of order 2 and order 3 polynomials with 2 points")
    {
        REQUIRE(
            integrate( [](double x) { return 3*x*x + 2*x + 1; },
                       legendre_rule<double, 2> ) ==
            doctest::Approx(4.0)
        );

        REQUIRE(
            integrate( [](double x) { return 4*x*x*x + 3*x*x + 2*x + 1; },
                       legendre_rule<double, 2> ) ==
            doctest::Approx(4.0)
        );
    }

    SUBCASE("Check the integral of an order 4 polynomial with 3 points")
    {
        constexpr auto f = [](double x)
        {
            return x*x*x*x + x*x*x + x*x + x + 1;
        };

        constexpr auto F = [](double x)
        {
            return x*x*x*x*x / 5 + x*x*x*x / 4 + x*x*x / 3 + x*x / 2 + x;
        };

        REQUIRE(integrate(f, legendre_rule<double, 3>) == doctest::Approx(F(1.0) - F(-1.0)));
    }
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

} /* namespace Quadrature */

} /* namespace Galerkin */

#endif /* QUADRATURE_HPP */
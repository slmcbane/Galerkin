/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef POLYNOMIALS_HPP
#define POLYNOMIALS_HPP

#include "FunctionBase.hpp"
#include "Metanomials.hpp"
#include "utils.hpp"

#include <array>

namespace Galerkin
{

namespace Polynomials
{

template <class T, class... Powers>
class Polynomial : public Functions::FunctionBase<Polynomial<T, Powers...>>
{
public:
    template <class... Args>
    constexpr Polynomial(Args... args) noexcept : m_coeffs{static_cast<T>(args)...}
    {
        static_assert(sizeof...(Args) == sizeof...(Powers));
    }

    constexpr Polynomial() : m_coeffs { 0 } {}

    constexpr const T& operator[](int i) const noexcept { return m_coeffs[i]; }
    T& operator[](int i) noexcept { return m_coeffs[i]; }

    template <class X>
    constexpr auto operator()(const X &args) const noexcept
    {
        return static_sum<0, sizeof...(Powers)>(
            [&](auto I) {
                return m_coeffs[I()] * Metanomials::raise(args, get<I()>(typeconst_list<Powers...>()));
            },
            zero<T>
        );
    }

    constexpr auto& coeffs() const noexcept { return m_coeffs; }

private:
    std::array<T, sizeof...(Powers)> m_coeffs;
};

template <class T, class... Powers>
constexpr bool operator==(const Polynomial<T, Powers...> &p1,
                          const Polynomial<T, Powers...> &p2) noexcept
{
    return p1.coeffs() == p2.coeffs();
}

template <class T, class... P1s, class... P2s>
constexpr bool operator==(const Polynomial<T, P1s...>&, const Polynomial<T, P2s...>&) noexcept
{
    return false;
}

/********************************************************************************
 * Test polynomial construction, arithmetic operations, evaluation.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Polynomials] Basic arithmetic and evaluation of polynomials")
{
    SUBCASE("A single-variable polynomial test case")
    {
        // p(x) = x^2 - 1
        constexpr auto p = Polynomial<int, Metanomials::Powers<2>, Metanomials::Powers<0>>(1, -1);

        // g(x) = x + 1;
        constexpr auto g = Polynomial<int, Metanomials::Powers<1>, Metanomials::Powers<0>>(1, 1);

        REQUIRE(p(std::tuple(1)) == 0);
        REQUIRE(p(std::tuple(1.5)) == doctest::Approx(1.25));
        REQUIRE(p(std::tuple(2)) == 3);

        REQUIRE(g(std::tuple(1)) == 2);
        REQUIRE(g(std::tuple(1.5)) == doctest::Approx(2.5));
        REQUIRE(g(std::tuple(2)) == 3);

        // Should work with an array argument, too.
        constexpr std::array x{3.0};
        REQUIRE(g(x) == doctest::Approx(4.0));

        // Check sum, product, and quotient of polynomials.
        constexpr auto p_plus_g = p + g;
        REQUIRE(p_plus_g(std::tuple(1)) == 2);
        REQUIRE(p_plus_g(std::tuple(1.5)) == doctest::Approx(3.75));
        REQUIRE(p_plus_g(std::tuple(2)) == 6);
        REQUIRE(p_plus_g(x) == doctest::Approx(12.0));

        constexpr auto p_times_g = p * g;
        REQUIRE(p_times_g(std::tuple(1)) == 0);
        REQUIRE(p_times_g(std::tuple(1.5)) == doctest::Approx(1.25 * 2.5));
        REQUIRE(p_times_g(std::tuple(2)) == 9);
        REQUIRE(p_times_g(x) == doctest::Approx(32.0));

        constexpr auto p_by_g = p / g;
        REQUIRE(p_by_g(std::tuple(1)) == 0);
        REQUIRE(p_by_g(std::tuple(1.5)) == doctest::Approx(1.25 / 2.5));
        REQUIRE(p_by_g(std::tuple(2)) == 1);
        REQUIRE(p_by_g(x) == doctest::Approx(2.0));
    }

    SUBCASE("Test functionality with a multi-variable polynomial")
    {
        // p(x, y) = x^2 + 2xy - 2y + 1
        constexpr auto p = Polynomial<double, Metanomials::Powers<2, 0>,
                                      Metanomials::Powers<1, 1>,
                                      Metanomials::Powers<0, 1>,
                                      Metanomials::Powers<0, 0>
                                     >(1, 2, -2, 1);
        // g(x, y) = y^2 - x + 2
        constexpr auto g = Polynomial<double, Metanomials::Powers<0, 2>,
                                      Metanomials::Powers<1, 0>, Metanomials::Powers<0, 0>
                                     >(1, -1, 2);

        // Evaluation points
        constexpr auto point1 = std::tuple(Rationals::rational<1, 2>, -Rationals::rational<1, 2>);
        constexpr auto point2 = std::array<double, 2>{3.0 / 4, 0 };
        constexpr auto point3 = std::tuple(1.0, -3.0 / 2);

        REQUIRE(p(point1) == doctest::Approx(7.0 / 4));
        REQUIRE(p(point2) == doctest::Approx(25.0 / 16));
        REQUIRE(p(point3) == doctest::Approx(2.0));

        REQUIRE(g(point1) == doctest::Approx(7.0 / 4));
        REQUIRE(g(point2) == doctest::Approx(5.0 / 4));
        REQUIRE(g(point3) == doctest::Approx(13.0 / 4));

        constexpr auto p_plus_g = p + g;
        constexpr auto p_times_g = p * g;
        constexpr auto p_by_g = p / g;

        REQUIRE(p_plus_g(point1) == doctest::Approx(14.0 / 4));
        REQUIRE(p_plus_g(point2) == doctest::Approx(25.0 / 16 + 5.0 / 4));
        REQUIRE(p_plus_g(point3) == doctest::Approx(2.0 + 13.0 / 4));

        REQUIRE(p_times_g(point1) == doctest::Approx(49.0 / 16));
        REQUIRE(p_times_g(point2) == doctest::Approx(125.0 / 64.0));
        REQUIRE(p_times_g(point3) == doctest::Approx(26.0 / 4));

        REQUIRE(p_by_g(point1) == doctest::Approx(1.0));
        REQUIRE(p_by_g(point2) == doctest::Approx(5.0 / 4));
        REQUIRE(p_by_g(point3) == doctest::Approx(8.0 / 13));
    }
}

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace Polynomials

} // namespace Galerkin

#endif /* POLYNOMIALS_HPP */
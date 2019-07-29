#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include "utils.hpp"

#include <array>
#include <cstdint>
#include <type_traits>

namespace Galerkin
{

namespace Legendre
{

constexpr auto gcd(uint64_t a, uint64_t b)
{
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

template <auto N>
constexpr auto gcd(std::array<uint64_t, N> numbers)
{
    if constexpr (N == 2)
    {
        return gcd(numbers[0], numbers[1]);
    }
    else
    {
        std::array<uint64_t, N - 1> newnumbers = {gcd(numbers[0], numbers[1])};
        for (int i = 2; i < N; ++i)
        {
            newnumbers[i - 1] = numbers[i];
        }
        return gcd(newnumbers);
    }
}

template <auto ORDER>
struct Polynomial
{
private:
    std::array<int64_t, ORDER + 1> m_coeffs;
    uint64_t m_denom;

public:
    typedef int64_t coeff_t;
    typedef uint64_t denom_t;

    constexpr Polynomial(std::array<coeff_t, ORDER + 1> c, denom_t d) noexcept : m_coeffs{0}, m_denom{1}
    {
        if (d != 0)
        {
            std::array<denom_t, ORDER + 2> ns = {0};
            for (int i = 0; i < ORDER + 1; ++i)
            {
                ns[i] = c[i] > 0 ? c[i] : -c[i];
            }
            ns[ORDER + 1] = d;
            auto div = gcd(ns);
            for (int i = 0; i < ORDER + 1; ++i)
            {
                // Note the static cast! You get promotion to the unsigned type
                // otherwise.
                m_coeffs[i] = c[i] / static_cast<coeff_t>(div);
            }
            m_denom = d / div;
        }
    }

    constexpr auto denominator() const noexcept { return m_denom; }
    constexpr auto coeffs() const noexcept { return m_coeffs; }

    template <class T>
    constexpr auto operator()(T x) const noexcept
    {
        T total = 0;
        for (int i = 0; i < ORDER + 1; ++i)
        {
            auto term = m_coeffs[i] * one<T>();
            for (int j = 0; j < i; ++j)
            {
                term = term * x;
            }
            total += term;
        }
        return total / m_denom;
    }

    constexpr auto derivative() const noexcept
    {
        if constexpr (ORDER == 0)
        {
            return Polynomial<0>({0}, 1);
        }
        else if constexpr (ORDER == 1)
        {
            return Polynomial<0>({m_coeffs[1]}, 1);
        }
        else
        {
            return static_sum<2, ORDER + 1>(
                [=](auto I) {
                    std::array<coeff_t, I> coeffs = {0};
                    coeffs[I - 1] = I * m_coeffs[I];
                    return Polynomial<I - 1>(coeffs, m_denom);
                },
                Polynomial<0>({m_coeffs[1]}, m_denom));
        }
    }
};

template <auto O1, auto O2>
constexpr auto operator+(Polynomial<O1> p1, Polynomial<O2> p2) noexcept
{
    const auto newdenom = (p1.denominator() * p2.denominator()) /
                          gcd(p1.denominator(), p2.denominator());

    if constexpr (O1 > O2)
    {
        auto coeffs = p1.coeffs();
        for (int i = 0; i < O2 + 1; ++i)
        {
            coeffs[i] = coeffs[i] * (newdenom / p1.denominator()) +
                        p2.coeffs()[i] * (newdenom / p2.denominator());
        }
        for (int i = O2 + 1; i < O1 + 1; ++i)
        {
            coeffs[i] = coeffs[i] * (newdenom / p1.denominator());
        }
        return Polynomial<O1>(coeffs, newdenom);
    }
    else
    {
        auto coeffs = p2.coeffs();
        for (int i = 0; i < O1 + 1; ++i)
        {
            coeffs[i] = coeffs[i] * (newdenom / p2.denominator()) +
                        p1.coeffs()[i] * (newdenom / p1.denominator());
        }
        for (int i = O1 + 1; i < O2 + 1; ++i)
        {
            coeffs[i] = coeffs[i] * (newdenom / p2.denominator());
        }
        return Polynomial<O2>(coeffs, newdenom);
    }
}

template <auto O1, auto O2>
constexpr bool operator==(Polynomial<O1> p1, Polynomial<O2> p2)
{
    if constexpr (O1 > O2)
    {
        std::array<int64_t, O1 + 1> coeffs{0};
        for (int i = 0; i < O2 + 1; ++i)
        {
            coeffs[i] = p2.coeffs()[i];
        }
        return p1 == Polynomial<O1>(coeffs, p2.denominator());
    }
    else if constexpr (O1 == O2)
    {
        return p1.coeffs() == p2.coeffs() && p1.denominator() == p2.denominator();
    }
    else
    {
        return p2 == p1;
    }
}

template <auto O1, auto O2>
constexpr auto operator*(Polynomial<O1> p1, Polynomial<O2> p2)
{
    constexpr auto new_order = O1 + O2;
    std::array<int64_t, new_order + 1> coeffs{0};
    for (int i = 0; i < O2 + 1; ++i)
    {
        for (int j = 0; j < O1 + 1; ++j)
        {
            coeffs[i + j] += p2.coeffs()[i] * p1.coeffs()[j];
        }
    }
    return Polynomial<new_order>(coeffs, p1.denominator() * p2.denominator());
}

template <auto O, class I>
constexpr std::enable_if_t<std::is_integral_v<I>, Polynomial<O>>
operator/(Polynomial<O> p, I v)
{
    return Polynomial<O>(p.coeffs(), v * p.denominator());
}

template <auto N>
constexpr auto polynomial =
    (Polynomial<1>({0, 2 * N - 1}, 1) * polynomial<N - 1> +
     Polynomial<0>({-(N - 1)}, 1) * polynomial<N - 2>) /
    N;

template <>
constexpr auto polynomial<0> = Polynomial<0>({1}, 1);

template <>
constexpr auto polynomial<1> = Polynomial<1>({0, 1}, 1);

template <class T, auto O>
constexpr T interval_root(Polynomial<O> p, T low, T high)
{
    constexpr auto pprime = p.derivative();
    auto x = (high + low) / 2;

    while (p(x) != zero<T>())
    {
        auto delta = -(p(x) / pprime(x));
        if (x + delta == x)
        {
            break;
        }
        x += delta;
    }
    return x;
}

template <class T, auto O>
constexpr auto all_roots(Polynomial<O> p)
{
    if constexpr (O == 1)
    {
        return std::array<T, 1> { -p(zero<T>()) / p.derivative()(zero<T>()) };
    }
    else
    {
        constexpr auto inflection_points = all_roots<T, O-1>(p.derivative());
        std::array<T, O> roots = {0};
        roots[0] = interval_root<T, O>(p, -one<T>(), inflection_points[0]);
        for (int i = 0; i < O - 2; ++i)
        {
            roots[i] = interval_root<T, O>(p, inflection_points[i], inflection_points[i + 1]);
        }
        roots[O - 1] = interval_root<T, O>(p, inflection_points[O - 2], one<T>());
        return roots;
    }
}

template <class T, auto O>
constexpr auto roots = all_roots<T, O>(polynomial<O>);

} /* namespace Legendre */

} /* namespace Galerkin */

/********************************************************************************
 * Start tests.
 * Please take the test code to be example of the API - this allows to keep
 * doc comments short.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using namespace Galerkin::Legendre;
using coeff_t = int64_t;

TEST_CASE("[Galerkin::Legendre] Test polynomial simplification")
{
    constexpr Polynomial<2> p1({1, 1, 1}, 2);
    constexpr Polynomial<3> p2({1, 0, 0, 1}, 2);
    constexpr Polynomial<2> p3({2, 0, 4}, 4);
    constexpr Polynomial<5> p4({6, 6, 9, 12, 18, 0}, 6);

    SUBCASE("These polynomials can't be reduced further")
    {
        REQUIRE(p1.coeffs() == std::array<coeff_t, 3>({1, 1, 1}));
        REQUIRE(p1.denominator() == 2);

        REQUIRE(p2.coeffs() == std::array<coeff_t, 4>({1, 0, 0, 1}));
        REQUIRE(p2.denominator() == 2);
    }

    SUBCASE("These polynomials should be reduced")
    {
        REQUIRE(p3 == Polynomial<2>({1, 0, 2}, 2));
        REQUIRE(p4 == Polynomial<5>({2, 2, 3, 4, 6, 0}, 2));
    }
}

TEST_CASE("[Galerkin::Legendre] Test polynomial addition")
{
    constexpr Polynomial<2> p1({1, 1, 1}, 2);
    constexpr Polynomial<3> p2({1, 0, 0, 1}, 2);
    constexpr Polynomial<2> p3({2, 0, 4}, 3);
    constexpr Polynomial<5> p4({6, 6, 9, 12, 18, 0}, 4);

    constexpr Polynomial<5> p5 = p1 + p2 + p3 + p4;
    REQUIRE(p5 == Polynomial<5>({38, 24, 49, 42, 54, 0}, 12));

    constexpr auto p6 = p5 + Polynomial<2>({1, 0, 2}, 12);
    REQUIRE(p6 == Polynomial<5>({13, 8, 17, 14, 18, 0}, 4));
}

TEST_CASE("[Galerkin::Legendre] Test polynomial multiplication")
{
    constexpr Polynomial<2> p1({2, 5, 4}, 3);
    constexpr Polynomial<3> p2({1, 2, 3, 4}, 4);

    constexpr auto p3 = p1 * p2;
    REQUIRE(p3 == Polynomial<5>({2, 9, 20, 31, 32, 16}, 12));
    REQUIRE(p3 == p2 * p1);
}

TEST_CASE("[Galerkin::Legendre] Test polynomial division by a scalar")
{
    constexpr Polynomial<3> p1({2, 6, 4, 2}, 7);
    REQUIRE(p1 / 2 == Polynomial<3>({1, 3, 2, 1}, 7));
}

TEST_CASE("[Galerkin::Legendre] Test polynomial differentiation")
{
    constexpr Polynomial<3> p({1, 2, 3, 2}, 2);
    REQUIRE(p.derivative() == Polynomial<2>({1, 3, 3}, 1));

    constexpr Polynomial<4> p2({1, 3, 5, 3, 1}, 3);
    REQUIRE(p2.derivative() == Polynomial<3>({3, 10, 9, 4}, 3));
}

TEST_CASE("[Galerkin::Legendre] Test correctness of Legendre::polynomial")
{
    REQUIRE(polynomial<2> == Polynomial<2>({-1, 0, 3}, 2));
    REQUIRE(polynomial<3> == Polynomial<3>({0, -3, 0, 5}, 2));
    REQUIRE(polynomial<7> ==
            Polynomial<7>({0, -35, 0, 315, 0, -693, 0, 429}, 16));
    REQUIRE(polynomial<10> ==
            Polynomial<10>({-63, 0, 3465, 0, -30030, 0, 90090, 0, -109395, 0, 46189}, 256));
}

#endif

#endif /* LEGENDRE_HPP */
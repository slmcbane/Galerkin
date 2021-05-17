/*
Copyright 2020 Sean McBane

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef POLYNOMIAL_POLYNOMIAL_HPP
#define POLYNOMIAL_POLYNOMIAL_HPP

#include "FunctionBase.hpp"
#include "Powers.hpp"

#include <array>
#include <type_traits>
#include <utility>

namespace Galerkin
{

namespace Polynomials
{

namespace detail
{

template <class T, class Tuple, std::size_t... Is>
constexpr std::array<T, sizeof...(Is)> to_array_impl(std::index_sequence<Is...>, const Tuple &tup) noexcept
{
    return std::array<T, sizeof...(Is)>{std::get<Is>(tup)...};
}

template <class T, class... Ts>
constexpr std::array<T, sizeof...(Ts)> to_array(const std::tuple<Ts...> &tup) noexcept
{
    return to_array_impl<T>(std::make_index_sequence<sizeof...(Ts)>(), tup);
}

template <class T, std::size_t N>
constexpr std::array<T, N> initialize_array()
{
    if constexpr (std::is_class_v<T>)
    {
        return std::array<T, N>();
    }
    else
    {
        return std::array<T, N>{0};
    }
}

struct PolyMaker;

template <std::size_t... Is, class... Xs, class... Ps, class T>
constexpr auto eval_impl(
    const std::array<T, sizeof...(Is)> &coeffs, std::index_sequence<Is...>, PowersList<Ps...>,
    const Xs &...xs) noexcept
{
    // using result_type = std::decay_t<decltype(std::get<0>(std::tuple(raise(Ps{}, xs...)...)))>;
    return ((raise(Ps{}, xs...) * coeffs[Is]) + ...);
}

} // namespace detail

template <class T, class... Ps>
class Polynomial : public Functions::FunctionBase<Polynomial<T, Ps...>>
{
    std::array<T, sizeof...(Ps)> m_coeffs;

    constexpr Polynomial(const std::array<T, sizeof...(Ps)> &cs) noexcept : m_coeffs(cs) {}

    template <class... Ts, class = std::enable_if_t<sizeof...(Ts) == sizeof...(Ps)>>
    constexpr Polynomial(const std::tuple<Ts...> &cs) noexcept : m_coeffs{detail::to_array<T>(cs)}
    {
    }

    friend struct detail::PolyMaker;

    template <std::size_t... Is, unsigned... As, class... Qs>
    constexpr auto partial_impl(
        std::index_sequence<Is...>, std::integer_sequence<unsigned, As...>,
        PowersList<Qs...>) const noexcept
    {
        if constexpr (sizeof...(Is) == 0)
        {
            return Polynomial();
        }
        else
        {
            const auto new_coeffs = std::array{(m_coeffs[Is] * As)...};
            return make_poly(new_coeffs, PowersList<Qs...>{});
        }
    }

  public:
    constexpr const auto &coeffs() const noexcept { return m_coeffs; }
    static constexpr auto num_terms = sizeof...(Ps);

    constexpr Polynomial() noexcept : m_coeffs{T(0)} {}

    constexpr bool operator==(const Polynomial<T, Ps...> &other) const noexcept
    {
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            if (!(other.m_coeffs[i] == m_coeffs[i]))
            {
                return false;
            }
        }
        return true;
    }

    constexpr Polynomial<T, Ps...> operator+(const Polynomial<T, Ps...> &other) const noexcept
    {
        Polynomial<T, Ps...> result;
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            result.m_coeffs[i] = m_coeffs[i] + other.m_coeffs[i];
        }
        return result;
    }

    template <class U>
    constexpr Polynomial<std::common_type_t<T, U>, Ps...> operator*(U x) const noexcept
    {
        auto new_coeffs = detail::initialize_array<std::common_type_t<T, U>, sizeof...(Ps)>();
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            new_coeffs[i] = m_coeffs[i] * x;
        }
        return make_poly(new_coeffs, PowersList<Ps...>{});
    }

    template <class F, class = std::enable_if_t<!std::is_same_v<F, Polynomial<T, Ps...>>>>
    constexpr auto operator*(const Functions::FunctionBase<F> &f) const
    {
        return Functions::FunctionBase<Polynomial<T, Ps...>>::operator*(static_cast<const F &>(f));
    }

    template <class F>
    constexpr auto operator/(const Functions::FunctionBase<F> &f) const
    {
        return Functions::FunctionBase<Polynomial<T, Ps...>>::operator/(static_cast<const F &>(f));
    } 

    template <class U>
    constexpr Polynomial<std::common_type_t<T, U>, Ps...> operator/(U x) const noexcept
    {
        auto new_coeffs = detail::initialize_array<std::common_type_t<T, U>, sizeof...(Ps)>();
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            new_coeffs[i] = m_coeffs[i] / x;
        }
        return make_poly(new_coeffs, PowersList<Ps...>{});
    }

    template <class U>
    constexpr Polynomial<T, Ps...> &operator*=(U x)
    {
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            m_coeffs[i] *= x;
        }
        return *this;
    }

    constexpr Polynomial<T, Ps...> &operator+=(const Polynomial<T, Ps...> &other) noexcept
    {
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            m_coeffs[i] += other.m_coeffs[i];
        }
        return *this;
    }

    template <class... Xs>
    constexpr auto operator()(const Xs &...xs) const noexcept
    {
        return detail::eval_impl(
            m_coeffs, std::make_index_sequence<num_terms>(), PowersList<Ps...>{}, xs...);
    }

    template <std::size_t I>
    constexpr auto partial() const noexcept
    {
        constexpr auto tup = partials_with_multipliers<I>(PowersList<Ps...>{});
        constexpr auto indices = std::get<0>(tup);
        constexpr auto constants = std::get<1>(tup);
        constexpr auto powers = std::get<2>(tup);
        return partial_impl(indices, constants, powers);
    }
};

namespace detail
{

struct PolyMaker
{
    template <class T, class... Ps>
    constexpr static auto create(const std::array<T, sizeof...(Ps)> &coeffs, PowersList<Ps...>) noexcept
    {
        return Polynomial<T, Ps...>(coeffs);
    }

    template <class... Ts, class... Ps>
    constexpr static auto create(const std::tuple<Ts...> &coeffs, PowersList<Ps...>) noexcept
    {
        static_assert(
            sizeof...(Ts) == sizeof...(Ps), "Length of coeffs should match number of Powers terms");
        using T = std::common_type_t<std::decay_t<Ts>...>;

        return Polynomial<T, Ps...>(coeffs);
    }
};

template <std::size_t FinalSize, class C, std::size_t... Is, std::size_t... Js>
constexpr auto
collect_coeffs_impl(const C &coeffs, std::index_sequence<Is...>, std::index_sequence<Js...>) noexcept
{
    using T = std::decay_t<decltype(coeffs[0])>;
    std::array<T, FinalSize> collected = initialize_array<T, FinalSize>();
    ((collected[Js] += coeffs[Is]), ...);
    return collected;
}

template <class C, std::size_t... Is, std::size_t FinalSize>
constexpr auto collect_coeffs(
    const C &coeffs, std::index_sequence<Is...> mapped_indices,
    std::integral_constant<std::size_t, FinalSize>) noexcept
{
    return collect_coeffs_impl<FinalSize>(
        coeffs, std::make_index_sequence<sizeof...(Is)>(), mapped_indices);
}

template <std::size_t FinalSize, class... Ts, std::size_t... Is, std::size_t... Js>
constexpr auto collect_coeffs_impl(
    const std::tuple<Ts...> &coeffs, std::index_sequence<Is...>, std::index_sequence<Js...>) noexcept
{
    auto collected = initialize_array<std::common_type_t<Ts...>, FinalSize>();
    ((collected[Js] += std::get<Is>(coeffs)), ...);
    return collected;
}

template <class... Ts, std::size_t... Is, std::size_t FinalSize>
constexpr auto collect_coeffs(
    const std::tuple<Ts...> &coeffs, std::index_sequence<Is...> mapped_indices,
    std::integral_constant<std::size_t, FinalSize>) noexcept
{
    static_assert(sizeof...(Ts) == sizeof...(Is));
    return collect_coeffs_impl<FinalSize>(
        coeffs, std::make_index_sequence<sizeof...(Ts)>(), mapped_indices);
}

template <class C, class... Ps>
struct size_checker : public std::false_type
{
};

template <class T, std::size_t N, class... Ps>
struct size_checker<T[N], Ps...> : public std::bool_constant<N == sizeof...(Ps)>
{
};

template <class T, std::size_t N, class... Ps>
struct size_checker<std::array<T, N>, Ps...> : public std::bool_constant<N == sizeof...(Ps)>
{
};

template <class... Ts, class... Ps>
struct size_checker<std::tuple<Ts...>, Ps...> : public std::bool_constant<sizeof...(Ts) == sizeof...(Ps)>
{
};

} // namespace detail

template <class C, class... Ps>
constexpr auto make_poly(const C &coeffs, PowersList<Ps...>) noexcept
{
    static_assert(detail::size_checker<C, Ps...>::value, "Wrong number of coefficients to make_poly");
    constexpr auto inds_and_powers = unique_and_sorted(PowersList<Ps...>{});
    constexpr auto mapped_indices = inds_and_powers.first;
    constexpr auto final_powers = inds_and_powers.second;

    return detail::PolyMaker::create(
        detail::collect_coeffs(
            coeffs, mapped_indices, std::integral_constant<std::size_t, final_powers.size>{}),
        final_powers);
}

template <class T, class... Ps, class U, class... Qs>
constexpr auto operator+(const Polynomial<T, Ps...> &p, const Polynomial<U, Qs...> &q) noexcept
{
    using V = std::common_type_t<T, U>;
    auto coeffs = detail::initialize_array<V, sizeof...(Ps) + sizeof...(Qs)>();
    for (unsigned i = 0; i < sizeof...(Ps); ++i)
    {
        coeffs[i] = p.coeffs()[i];
    }
    for (unsigned i = 0; i < sizeof...(Qs); ++i)
    {
        coeffs[i + p.num_terms] = q.coeffs()[i];
    }
    return make_poly(coeffs, PowersList<Ps..., Qs...>{});
}

template <class T, class... Ps, class U, class... Qs>
constexpr auto operator*(const Polynomial<T, Ps...> &p, const Polynomial<U, Qs...> &q) noexcept
{
    using V = decltype(std::declval<T>() * std::declval<U>());
    auto coeffs = detail::initialize_array<V, sizeof...(Ps) * sizeof...(Qs)>();
    for (unsigned i = 0; i < sizeof...(Ps); ++i)
    {
        for (unsigned j = 0; j < sizeof...(Qs); ++j)
        {
            coeffs[i * sizeof...(Qs) + j] = p.coeffs()[i] * q.coeffs()[j];
        }
    }

    return make_poly(coeffs, PowersList<Ps...>{} * PowersList<Qs...>{});
}

template <class T, class U, class... Ps>
constexpr auto operator*(U x, const Polynomial<T, Ps...> &p) noexcept
{
    return p * x;
}

template <std::size_t I, class T, class... Ps>
constexpr auto partial(const Polynomial<T, Ps...> &p) noexcept
{
    return p.template partial<I>();
}

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Add two single-variable polynomials with same powers")
{
    constexpr auto coeffs1 = std::tuple(1, 2, 3, 4);
    constexpr auto coeffs2 = std::tuple(5, 6, 7, 8);
    constexpr auto powers = PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{};
    constexpr auto poly1 = make_poly(coeffs1, powers);
    constexpr auto poly2 = make_poly(coeffs2, powers);
    constexpr auto result = poly1 + poly2;

    REQUIRE(result.num_terms == 4);
    REQUIRE(result.coeffs()[0] == 6);
    REQUIRE(result.coeffs()[1] == 8);
    REQUIRE(result.coeffs()[2] == 10);
    REQUIRE(result.coeffs()[3] == 12);
}

TEST_CASE("Add two single-variable polynomials with same powers, in place")
{
    constexpr auto coeffs1 = std::tuple(1, 2, 3, 4);
    constexpr auto coeffs2 = std::tuple(5, 6, 7, 8);
    constexpr auto powers = PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{};
    auto poly1 = make_poly(coeffs1, powers);
    constexpr auto poly2 = make_poly(coeffs2, powers);
    poly1 += poly2;

    REQUIRE(poly1.coeffs()[0] == 6);
    REQUIRE(poly1.coeffs()[1] == 8);
    REQUIRE(poly1.coeffs()[2] == 10);
    REQUIRE(poly1.coeffs()[3] == 12);
}

TEST_CASE("Add two single-variable polynomials with different powers")
{
    constexpr auto coeffs1 = std::tuple(1, 2, 3, 4);
    constexpr auto coeffs2 = std::tuple(5, 6, 7, 8);
    constexpr auto powers1 = PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{};
    constexpr auto powers2 = PowersList<Powers<2>, Powers<3>, Powers<4>, Powers<5>>{};
    constexpr auto poly1 = make_poly(coeffs1, powers1);
    constexpr auto poly2 = make_poly(coeffs2, powers2);
    constexpr auto result = poly1 + poly2;

    REQUIRE(result.num_terms == 6);
    REQUIRE(result.coeffs()[0] == 1);
    REQUIRE(result.coeffs()[1] == 2);
    REQUIRE(result.coeffs()[2] == 8);
    REQUIRE(result.coeffs()[3] == 10);
    REQUIRE(result.coeffs()[4] == 7);
    REQUIRE(result.coeffs()[5] == 8);
}

TEST_CASE("Add two muli-variable polynomials with different powers")
{
    constexpr auto coeffs1 = std::tuple(2, -3, 1);
    constexpr auto coeffs2 = std::tuple(3, 4);
    constexpr auto powers1 = PowersList<Powers<1, 0>, Powers<1, 1>, Powers<0, 1>>{};
    constexpr auto powers2 = PowersList<Powers<0, 1>, Powers<2, 0>>{};

    constexpr auto poly1 = make_poly(coeffs1, powers1);
    constexpr auto poly2 = make_poly(coeffs2, powers2);

    constexpr auto result = poly1 + poly2;
    REQUIRE(result.num_terms == 4);
    REQUIRE(result.coeffs()[0] == 4);
    REQUIRE(result.coeffs()[1] == 2);
    REQUIRE(result.coeffs()[2] == -3);
    REQUIRE(result.coeffs()[3] == 4);
}

TEST_CASE("Construct a polynomial with coefficient array")
{
    SUBCASE("Only a single coefficient")
    {
        constexpr auto powers = PowersList<Powers<0>>{};
        constexpr auto coeffs = std::array<double, 1>{1.0};
        constexpr auto poly = make_poly(coeffs, powers);
        REQUIRE(poly.coeffs()[0] == 1.0);
    }

    SUBCASE("Single-variable polynomial, powers are sorted, no duplicates")
    {
        constexpr auto powers = PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{};
        constexpr std::array<double, 4> coeffs{1.0, 2.0, 3.0, 4.0};
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(poly.coeffs()[0] == 1.0);
        REQUIRE(poly.coeffs()[1] == 2.0);
        REQUIRE(poly.coeffs()[2] == 3.0);
        REQUIRE(poly.coeffs()[3] == 4.0);
    }

    SUBCASE("Single-variable polynomial, powers are not sorted, no duplicates")
    {
        constexpr auto powers = PowersList<Powers<0>, Powers<3>, Powers<2>, Powers<1>>{};
        constexpr std::array<double, 4> coeffs{1.0, 2.0, 3.0, 4.0};
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(poly.coeffs()[0] == 1.0);
        REQUIRE(poly.coeffs()[1] == 4.0);
        REQUIRE(poly.coeffs()[2] == 3.0);
        REQUIRE(poly.coeffs()[3] == 2.0);
    }

    SUBCASE("Single-variable polynomial, unsorted and duplicate powers")
    {
        constexpr auto powers = PowersList<
            Powers<3>, Powers<4>, Powers<4>, Powers<2>, Powers<1>, Powers<3>, Powers<2>, Powers<0>,
            Powers<1>>{};
        constexpr std::array<int, powers.size> coeffs{3, 4, 5, 4, 3, 5, 3, 5, 3};
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(poly.num_terms == 5);

        REQUIRE(poly.coeffs()[0] == 5);
        REQUIRE(poly.coeffs()[1] == 6);
        REQUIRE(poly.coeffs()[2] == 7);
        REQUIRE(poly.coeffs()[3] == 8);
        REQUIRE(poly.coeffs()[4] == 9);
    }

    SUBCASE("Multiple variables, unsorted and duplicate powers")
    {
        constexpr auto powers = PowersList<
            Powers<1, 0>, Powers<0, 0>, Powers<1, 0>, Powers<0, 1>, Powers<1, 1>, Powers<0, 1>>{};
        constexpr std::array<int, powers.size> coeffs{2, 3, 4, 5, 6, 7};
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(poly.num_terms == 4);

        REQUIRE(poly.coeffs()[0] == 3);
        REQUIRE(poly.coeffs()[1] == 12);
        REQUIRE(poly.coeffs()[2] == 6);
        REQUIRE(poly.coeffs()[3] == 6);
    }
}

TEST_CASE("Construct a polynomial with a coefficient tuple")
{
    SUBCASE("Single-variable polynomial, powers are sorted, no duplicates")
    {
        constexpr auto powers = PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{};
        constexpr auto coeffs = std::tuple(1, 2, 3, '\4');
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(std::is_same_v<
                std::remove_cv_t<decltype(poly)>,
                Polynomials::Polynomial<int, Powers<0>, Powers<1>, Powers<2>, Powers<3>>>);

        REQUIRE(poly.coeffs()[0] == 1);
        REQUIRE(poly.coeffs()[1] == 2);
        REQUIRE(poly.coeffs()[2] == 3);
        REQUIRE(poly.coeffs()[3] == 4);
    }

    SUBCASE("Single-variable polynomial, powers are not sorted, no duplicates")
    {
        constexpr auto powers = PowersList<Powers<0>, Powers<3>, Powers<2>, Powers<1>>{};
        constexpr auto coeffs = std::tuple(1.0, 2.0, 3.0, 4.0);
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(poly.coeffs()[0] == 1.0);
        REQUIRE(poly.coeffs()[1] == 4.0);
        REQUIRE(poly.coeffs()[2] == 3.0);
        REQUIRE(poly.coeffs()[3] == 2.0);
    }

    SUBCASE("Single-variable polynomial, unsorted and duplicate powers")
    {
        constexpr auto powers = PowersList<
            Powers<3>, Powers<4>, Powers<4>, Powers<2>, Powers<1>, Powers<3>, Powers<2>, Powers<0>,
            Powers<1>>{};
        constexpr auto coeffs = std::tuple(3, 4, 5, 4, 3, 5, 3, 5, 3);
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(poly.num_terms == 5);

        REQUIRE(poly.coeffs()[0] == 5);
        REQUIRE(poly.coeffs()[1] == 6);
        REQUIRE(poly.coeffs()[2] == 7);
        REQUIRE(poly.coeffs()[3] == 8);
        REQUIRE(poly.coeffs()[4] == 9);
    }

    SUBCASE("Multiple variables, unsorted and duplicate powers")
    {
        constexpr auto powers = PowersList<
            Powers<1, 0>, Powers<0, 0>, Powers<1, 0>, Powers<0, 1>, Powers<1, 1>, Powers<0, 1>>{};
        constexpr auto coeffs = std::tuple(2, 3, 4, 5, 6, 7);
        constexpr auto poly = make_poly(coeffs, powers);

        REQUIRE(poly.num_terms == 4);

        REQUIRE(poly.coeffs()[0] == 3);
        REQUIRE(poly.coeffs()[1] == 12);
        REQUIRE(poly.coeffs()[2] == 6);
        REQUIRE(poly.coeffs()[3] == 6);
    }
}

TEST_CASE("Evaluating a multivariable polynomial")
{
    constexpr auto powers = PowersList<Powers<3, 0>, Powers<0, 3>, Powers<2, 1>, Powers<1, 2>,
        Powers<2, 0>, Powers<0, 2>, Powers<1, 1>, Powers<1, 0>, Powers<0, 1>, Powers<0, 0>>{};

    constexpr std::array<int, powers.size> coeffs { 2, 1, -3, 4, 0, 5, 1, 1, 2, 1 };

    constexpr auto poly = make_poly(coeffs, powers);

    REQUIRE(poly(1, 2) == 48);
    REQUIRE(poly(2, 2) == 63);
    REQUIRE(poly(4.0, '\2') == 141);
}

TEST_CASE("Squaring a single-variable polynomial")
{
    constexpr auto powers = PowersList<Powers<2>, Powers<1>, Powers<0>>{};
    constexpr auto coeffs = std::tuple(1, -1, 2);
    constexpr auto poly = make_poly(coeffs, powers);
    constexpr auto poly2 = poly * poly;

    static_assert(std::is_same_v<
        std::decay_t<decltype(poly2)>,
        Polynomials::Polynomial<int, Powers<0>, Powers<1>, Powers<2>, Powers<3>, Powers<4>>>);

    REQUIRE(poly2.coeffs()[0] == 4);
    REQUIRE(poly2.coeffs()[1] == -4);
    REQUIRE(poly2.coeffs()[2] == 5);
    REQUIRE(poly2.coeffs()[3] == -2);
    REQUIRE(poly2.coeffs()[4] == 1);
}

TEST_CASE("A more complicated, multi-variable test case")
{
    constexpr auto powers = PowersList<Powers<0, 0>, Powers<0, 1>, Powers<1, 0>, Powers<1, 1>>{};
    constexpr auto coeffs1 = std::array<int, 4>{ 2, 1, 2, 3 };
    constexpr auto coeffs2 = std::array<int, 4>{ 4, 3, 2, 1 };
    constexpr auto poly1 = make_poly(coeffs1, powers);
    constexpr auto poly2 = make_poly(coeffs2, powers);
    constexpr auto result = poly1 * poly2;

    static_assert(result.num_terms == 9);
    REQUIRE(result.coeffs()[0] == 8);
    REQUIRE(result.coeffs()[1] == 10);
    REQUIRE(result.coeffs()[2] == 3);
    REQUIRE(result.coeffs()[3] == 12);
    REQUIRE(result.coeffs()[4] == 22);
    REQUIRE(result.coeffs()[5] == 10);
    REQUIRE(result.coeffs()[6] == 4);
    REQUIRE(result.coeffs()[7] == 8);
    REQUIRE(result.coeffs()[8] == 3);
}

TEST_CASE("Scalar multiplication test")
{
    constexpr auto powers = PowersList<Powers<0, 0>, Powers<0, 1>, Powers<1, 0>, Powers<1, 1>>{};
    constexpr auto coeffs1 = std::array<int, 4>{2, 1, 2, 3};
    constexpr auto poly1 = make_poly(coeffs1, powers);
    auto poly2 = poly1 * 1.5;

    REQUIRE(poly2.coeffs()[0] == 3);
    REQUIRE(poly2.coeffs()[1] == 1.5);
    REQUIRE(poly2.coeffs()[2] == 3);
    REQUIRE(poly2.coeffs()[3] == 4.5);

    auto poly3 = poly1;
    poly3 *= 1.5;
    for (int i = 0; i < 4; ++i)
    {
        REQUIRE(poly3.coeffs()[i] == static_cast<int>(1.5) * poly3.coeffs()[i]);
    }
}

TEST_CASE("Partial derivatives of a single-variable function")
{
    constexpr auto powers = PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{};
    constexpr auto coefficients = std::array{1.0, 2.0, 3.0, 4.0};
    constexpr auto poly = make_poly(coefficients, powers);
    constexpr auto polyprime = partial<0>(poly);

    static_assert(
        std::is_same_v<
            std::remove_cv_t<decltype(polyprime)>, Polynomial<double, Powers<0>, Powers<1>, Powers<2>>>);

    REQUIRE(polyprime.coeffs()[0] == 2.0);
    REQUIRE(polyprime.coeffs()[1] == 6.0);
    REQUIRE(polyprime.coeffs()[2] == 12.0);
}

TEST_CASE("Partial derivatives for a function of two variables")
{
    constexpr auto powers = PowersList<
        Powers<0, 0>, Powers<0, 1>, Powers<0, 2>, Powers<1, 0>, Powers<1, 1>, Powers<1, 2>, Powers<2, 2>>{};
    constexpr int coefficients[] = {1, 2, 1, 2, 1, 2, 1};
    constexpr auto poly = make_poly(coefficients, powers);
    constexpr auto d0poly = partial<0>(poly);

    using T = decltype(std::declval<int>() * std::declval<unsigned>());
    static_assert(std::is_same_v<
        std::remove_cv_t<decltype(d0poly)>,
        Polynomial<T, Powers<0, 0>, Powers<0, 1>, Powers<0, 2>, Powers<1, 2>>>);

    REQUIRE(d0poly.coeffs()[0] == 2);
    REQUIRE(d0poly.coeffs()[1] == 1);
    REQUIRE(d0poly.coeffs()[2] == 2);
    REQUIRE(d0poly.coeffs()[3] == 2);

    constexpr auto d1poly = partial<1>(poly);
    static_assert(std::is_same_v<
        std::remove_cv_t<decltype(d1poly)>,
        Polynomial<T, Powers<0, 0>, Powers<0, 1>, Powers<1, 0>, Powers<1, 1>, Powers<2, 1>>>);
    REQUIRE(d1poly.coeffs()[0] == 2);
    REQUIRE(d1poly.coeffs()[1] == 2);
    REQUIRE(d1poly.coeffs()[2] == 1);
    REQUIRE(d1poly.coeffs()[3] == 4);
    REQUIRE(d1poly.coeffs()[4] == 2);
}

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace Polynomials

} // namespace Galerkin

#endif // POLYNOMIAL_POLYNOMIAL_HPP

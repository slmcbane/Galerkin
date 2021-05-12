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

struct PolyMaker;

template <std::size_t... Is, class... Xs, class... Ps, class T>
constexpr T eval_impl(
    const std::array<T, sizeof...(Is)> &coeffs, std::index_sequence<Is...>, PowersList<Ps...>,
    const Xs &...xs) noexcept
{
    return ((raise(Ps{}, xs...) * coeffs[Is]) + ...);
}

} // namespace detail

template <class T, class... Ps>
class Polynomial
{
    static_assert(std::is_arithmetic_v<T>);
    std::array<T, sizeof...(Ps)> m_coeffs;

    constexpr Polynomial(const std::array<T, sizeof...(Ps)> &cs) noexcept : m_coeffs(cs) {}

    constexpr Polynomial() noexcept : m_coeffs{0} {}

    template <class... Ts, class = std::enable_if_t<sizeof...(Ts) == sizeof...(Ps)>>
    constexpr Polynomial(const std::tuple<Ts...> &cs) noexcept : m_coeffs{detail::to_array<T>(cs)}
    {
    }

    friend struct detail::PolyMaker;

    template <std::size_t... Is, unsigned... As, class... Qs>
    constexpr auto partial_impl(
        std::index_sequence<Is...>, std::integer_sequence<unsigned, As...>, PowersList<Qs...>) const noexcept
    {
        const auto new_coeffs = std::array{(m_coeffs[Is] * As)...};
        return make_poly(new_coeffs, PowersList<Qs...>{});
    }

  public:
    constexpr const auto &coeffs() const noexcept { return m_coeffs; }
    static constexpr auto num_terms = sizeof...(Ps);

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
    constexpr std::enable_if_t<std::is_arithmetic_v<U>, Polynomial<std::common_type_t<T, U>, Ps...>>
    operator*(U x) const noexcept
    {
        std::array<std::common_type_t<T, U>, sizeof...(Ps)> new_coeffs{0};
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            new_coeffs[i] = m_coeffs[i] * x;
        }
        return make_poly(new_coeffs, PowersList<Ps...>{});
    }

    template <class U>
    constexpr std::enable_if_t<std::is_arithmetic_v<U>, Polynomial<std::common_type_t<T, U>, Ps...>>
    operator/(U x) const noexcept
    {
        std::array<std::common_type_t<T, U>, sizeof...(Ps)> new_coeffs{0};
        for (unsigned i = 0; i < sizeof...(Ps); ++i)
        {
            new_coeffs[i] = m_coeffs[i] / x;
        }
        return make_poly(new_coeffs, PowersList<Ps...>{});
    }

    template <class U>
    constexpr std::enable_if_t<std::is_arithmetic_v<U>, Polynomial<T, Ps...>> &operator*=(U x)
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
    constexpr T operator()(const Xs &...xs) const noexcept
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
    std::array<T, FinalSize> collected{0};
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
    std::array<std::common_type_t<Ts...>, FinalSize> collected{0};
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
    std::array<V, sizeof...(Ps) + sizeof...(Qs)> coeffs{0};
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
    std::array<V, sizeof...(Ps) * sizeof...(Qs)> coeffs{0};
    for (unsigned i = 0; i < sizeof...(Ps); ++i)
    {
        for (unsigned j = 0; j < sizeof...(Qs); ++j)
        {
            coeffs[i * sizeof...(Qs) + j] = p.coeffs()[i] * q.coeffs()[j];
        }
    }

    return make_poly(coeffs, PowersList<Ps...>{} * PowersList<Qs...>{});
}

template <std::size_t I, class T, class... Ps>
constexpr auto partial(const Polynomial<T, Ps...> &p) noexcept
{
    return p.template partial<I>();
}

} // namespace Polynomials

} // namespace Galerkin

#endif // POLYNOMIAL_POLYNOMIAL_HPP

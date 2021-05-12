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

#ifndef POLYNOMIAL_POWERS_HPP
#define POLYNOMIAL_POWERS_HPP

#include <array>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Galerkin
{

namespace Polynomials
{

template <unsigned... Ps>
struct Powers
{
    constexpr Powers() noexcept {}
    constexpr static auto sum = (Ps + ...);

    static constexpr std::array<unsigned, sizeof...(Ps)> terms{Ps...};

    template <int which>
    constexpr static unsigned term = terms[which];

    static constexpr std::size_t nvars = sizeof...(Ps);

    template <unsigned... Qs>
    constexpr auto operator*(Powers<Qs...>) const noexcept
    {
        static_assert(sizeof...(Ps) == sizeof...(Qs));
        return Powers<(Ps + Qs)...>{};
    }
};

template <unsigned... Ps, unsigned... Qs>
constexpr bool operator==(Powers<Ps...>, Powers<Qs...>) noexcept
{
    return std::is_same_v<Powers<Ps...>, Powers<Qs...>>;
}

template <unsigned... P1s, unsigned... P2s>
constexpr bool operator<(Powers<P1s...>, Powers<P2s...>) noexcept
{
    static_assert(
        sizeof...(P1s) == sizeof...(P2s),
        "Comparison of Powers types only meaningful if they are the same size");
    using Ps1 = Powers<P1s...>;
    using Ps2 = Powers<P2s...>;
    if constexpr (Ps1::sum < Ps2::sum)
    {
        return true;
    }
    else if constexpr (Ps1::sum > Ps2::sum)
    {
        return false;
    }
    else
    {
        constexpr unsigned terms1[sizeof...(P1s)] = {P1s...};
        constexpr unsigned terms2[sizeof...(P1s)] = {P2s...};
        for (unsigned i = 0; i < sizeof...(P1s); ++i)
        {
            if (terms1[i] < terms2[i])
            {
                return true;
            }
            else if (terms1[i] > terms2[i])
            {
                return false;
            }
        }
        return false;
    }
}

namespace detail
{

template <class T>
struct is_powers
{
    constexpr static bool value = false;
};

template <unsigned... Ps>
struct is_powers<Powers<Ps...>>
{
    constexpr static bool value = true;
};

template <class... Ps>
struct all_are_powers;

template <>
struct all_are_powers<>
{
    static constexpr bool value = true;
};

template <class P, class... Ps>
struct all_are_powers<P, Ps...>
{
    static constexpr bool value = is_powers<P>::value && all_are_powers<Ps...>::value;
};

template <class... Ps>
struct all_same_size;

template <>
struct all_same_size<> : public std::true_type
{
};

template <class P>
struct all_same_size<P> : public std::true_type
{
};

template <class P, class Q, class... Ps>
struct all_same_size<P, Q, Ps...>
    : public std::bool_constant<P::nvars == Q::nvars && all_same_size<Ps...>::value>
{
};

template <std::size_t... Is, class... Ps>
constexpr auto tuple_tail(std::index_sequence<Is...>, std::tuple<Ps...>)
{
    return std::tuple{std::get<Is + 2>(std::tuple{Ps{}...})...};
}

template <class... Ps>
constexpr auto tuple_tail(std::tuple<Ps...>)
{
    constexpr auto Is = std::make_index_sequence<sizeof...(Ps) - 2>();
    return tuple_tail(Is, std::tuple{Ps{}...});
}

} // namespace detail

template <class... Ps>
struct PowersList
{
    static_assert(detail::all_are_powers<Ps...>::value, "PowersList can only contain Powers objects");
    static_assert(
        detail::all_same_size<Ps...>::value,
        "All Powers in PowersList should have same order (number of variables)");
    constexpr static auto size = sizeof...(Ps);
    static constexpr std::size_t nvars = std::tuple_element_t<0, std::tuple<Ps...>>::nvars;

    template <class... Qs>
    constexpr auto operator+(PowersList<Qs...>) const noexcept
    {
        return PowersList<Ps..., Qs...>{};
    }

    template <std::size_t which>
    static constexpr auto term()
    {
        return std::tuple_element_t<which, std::tuple<Ps...>>{};
    }
};

template <unsigned P, class T>
constexpr T raise([[maybe_unused]] T x) noexcept
{
    if constexpr (P == 0)
    {
        return 1;
    }
    else if constexpr (P == 1)
    {
        return x;
    }
    else
    {
        constexpr unsigned first_power = P / 2;
        constexpr unsigned second_power = P - first_power;
        return raise<first_power>(x) * raise<second_power>(x);
    }
}

namespace detail
{

template <class Tup, unsigned... Ps, std::size_t... Is>
constexpr auto raise_tuple_impl(const Tup &xs, Powers<Ps...>, std::index_sequence<Is...>) noexcept
{
    return (raise<Ps>(std::get<Is>(xs)) * ...);
}

} // namespace detail

template <class T, unsigned... Ps>
constexpr auto raise(Powers<Ps...>, const std::array<T, sizeof...(Ps)> &xs) noexcept
{
    return detail::raise_tuple_impl(xs, Powers<Ps...>{}, std::make_index_sequence<sizeof...(Ps)>());
}

template <class... Ts, unsigned... Ps>
constexpr auto raise(Powers<Ps...>, const std::tuple<Ts...> &xs) noexcept
{
    return detail::raise_tuple_impl(xs, Powers<Ps...>{}, std::make_index_sequence<sizeof...(Ps)>());
}

template <class... Ts, unsigned... Ps>
constexpr auto raise(Powers<Ps...>, Ts... xs) noexcept
{
    return (raise<Ps>(xs) * ...);
}

template <unsigned... Ps, class... Qs>
constexpr auto operator*(Powers<Ps...>, PowersList<Qs...>) noexcept
{
    return PowersList<decltype(Powers<Ps...>{} * Qs{})...>{};
}

template <class... Ps, class... Qs>
constexpr auto operator*(PowersList<Ps...>, PowersList<Qs...>) noexcept
{
    return ((Ps{} * PowersList<Qs...>{}) + ...);
}

namespace detail
{

template <unsigned... Ps>
constexpr auto to_array(Powers<Ps...>) noexcept
{
    return std::array<unsigned, sizeof...(Ps)>{Ps...};
}

template <class P, class... Ps>
constexpr auto expand_powers(P, Ps...) noexcept
{
    return std::array<std::decay_t<decltype(to_array(P{}))>, sizeof...(Ps)+1>{to_array(P{}), to_array(Ps{})...};
}

template <std::size_t N, std::size_t... Is>
constexpr bool compare_arrays_impl(
    const std::array<unsigned, N> &a, const std::array<unsigned, N> &b, std::index_sequence<Is...>) noexcept
{
    return ((a[Is] == b[Is]) && ...);
}

template <std::size_t N>
constexpr bool compare_arrays(const std::array<unsigned, N> &a, const std::array<unsigned, N> &b) noexcept
{
    return compare_arrays_impl(a, b, std::make_index_sequence<N>());
}

template <std::size_t N>
constexpr bool array_less_than(const std::array<unsigned, N> &a, const std::array<unsigned, N> &b) noexcept
{
    for (std::size_t i = 0; i < N; ++i)
    {
        if (a[i] > b[i])
        {
            return false;
        }
        else if (a[i] < b[i])
        {
            return true;
        }
    }
    return false;
}

template <std::size_t n, std::size_t sz>
constexpr auto sort_array(const std::array<std::array<unsigned, n>, sz> &a) noexcept
{
    std::array<std::array<unsigned, n>, sz> cp = a;

    constexpr auto sorting_pass = [](auto &arr) {
        if (arr.size() < 2)
        {
            return false;
        }

        bool modified = false;
        for (auto it = arr.begin(); it != arr.end() - 1; ++it)
        {
            if (array_less_than(*(it + 1), *it))
            {
                auto tmp = *it;
                *it = *(it + 1);
                *(it + 1) = tmp;
                modified = true;
            }
        }
        return modified;
    };

    while (sorting_pass(cp))
        ;
    return cp;
}

template <class P, class A>
constexpr std::size_t map_power(const A &arr) noexcept
{
    constexpr auto to_compare = to_array(P{});
    for (std::size_t i = 0; i < arr.size(); ++i)
    {
        if (compare_arrays(to_compare, arr[i]))
        {
            return i;
        }
    }
}

template <std::size_t N, std::size_t Sz>
constexpr std::size_t num_unique_powers(const std::array<std::array<unsigned, N>, Sz> &arr) noexcept
{
    // Input array is already sorted.
    std::size_t count = 0;
    std::size_t index = 0;
    while (index < Sz)
    {
        count += 1;
        const auto val = arr[index];
        do
        {
            index += 1;
        } while (index < Sz && compare_arrays(arr[index], val));
    }
    return count;
}

template <std::size_t FinalSize, std::size_t N, std::size_t Sz>
constexpr std::array<std::array<unsigned, N>, FinalSize>
get_unique_array(const std::array<std::array<unsigned, N>, Sz> &arr) noexcept
{
    std::array<std::array<unsigned, N>, FinalSize> result{0};
    std::size_t index = 0;
    std::size_t count = 0;
    while (index < Sz)
    {
        result[count++] = arr[index];
        const auto val = arr[index];
        do
        {
            index += 1;
        } while (index < Sz && compare_arrays(arr[index], val));
    }
    return result;
}

template <class... Ps>
struct UniqueAndSorted
{
    static constexpr auto sorted_array = sort_array(expand_powers(Ps{}...));
    static constexpr auto value = get_unique_array<num_unique_powers(sorted_array)>(sorted_array);
    static constexpr auto size = value.size();
};

template <std::size_t i, class... Ps>
constexpr std::size_t get_power_index(PowersList<Ps...>) noexcept
{
    static_assert(i < UniqueAndSorted<Ps...>::size);
    constexpr auto expanded = expand_powers(Ps{}...);
    for (std::size_t j = 0; j < sizeof...(Ps); ++j)
    {
        if (compare_arrays(UniqueAndSorted<Ps...>::value[i], expanded[j]))
        {
            return j;
        }
    }
}

template <std::size_t... Is, class... Ps>
constexpr auto get_final_powers_impl(std::index_sequence<Is...>, PowersList<Ps...>) noexcept
{
    return PowersList<std::decay_t<decltype(PowersList<Ps...>::template term<get_power_index<Is>(PowersList<Ps...>{})>())>...>{};
}

template <class... Ps>
constexpr auto get_final_powers(PowersList<Ps...>) noexcept
{
    constexpr std::size_t num_unique = UniqueAndSorted<Ps...>::value.size();
    constexpr auto helper_seq = std::make_index_sequence<num_unique>();
    return get_final_powers_impl(helper_seq, PowersList<Ps...>{});
}

} // namespace detail

template <class... Ps>
constexpr auto unique_and_sorted(PowersList<Ps...>)
{
    constexpr auto final_indices =
        std::index_sequence<detail::map_power<Ps>(detail::UniqueAndSorted<Ps...>::value)...>{};
    constexpr auto final_powers = detail::get_final_powers(PowersList<Ps...>{});
    return std::make_pair(final_indices, final_powers);
}

namespace detail
{

template <std::size_t... Is, std::size_t... Js>
constexpr std::index_sequence<Is..., Js...> concatenate(std::index_sequence<Is...>, std::index_sequence<Js...>) noexcept
{
    return std::index_sequence<Is..., Js...>{};
}

template <std::size_t... Is>
constexpr std::index_sequence<(Is+1)...> add1to(std::index_sequence<Is...>) noexcept
{
    return std::index_sequence<(Is+1)...>{};
}

template <std::size_t I, class... Ps>
struct FilterZeroPowers;

template <std::size_t I>
struct FilterZeroPowers<I>
{
    typedef std::index_sequence<> seq_type;
    typedef PowersList<> plist_type;
};

template <std::size_t I, class P, class... Ps>
struct FilterZeroPowers<I, P, Ps...>
{
    typedef typename FilterZeroPowers<I, Ps...>::seq_type subseq_type;
    typedef typename FilterZeroPowers<I, Ps...>::plist_type sublist_type;

    constexpr static auto compute() noexcept
    {
        if constexpr (to_array(P{})[I] == 0)
        {
            return std::make_pair(add1to(subseq_type{}), sublist_type{});
        }
        else
        {
            return std::make_pair(concatenate(std::index_sequence<0>{}, add1to(subseq_type{})),
                PowersList<P>{} + sublist_type{});
        }
    }

    typedef decltype(compute().first) seq_type;
    typedef decltype(compute().second) plist_type;
};

template <std::size_t I, class... Ps>
constexpr auto get_partial_constants() noexcept
{
    return std::integer_sequence<unsigned, to_array(Ps{})[I]...>{};
}

template <std::size_t I, std::size_t... Js, class Ps>
constexpr auto apply_partial_impl(std::index_sequence<Js...>, Ps)
{
    return Powers<(Ps::terms[Js] - (Js == I))...>{};
}

template <std::size_t I, class Ps>
constexpr auto apply_partial(Ps) noexcept
{
    constexpr auto helper_sequence = std::make_index_sequence<Ps::nvars>();
    return apply_partial_impl<I>(helper_sequence, Ps{});
}

template <std::size_t I, class... Ps>
constexpr auto do_partials(PowersList<Ps...>) noexcept
{
    constexpr auto constants = get_partial_constants<I, Ps...>();
    constexpr PowersList<decltype(apply_partial<I>(Ps{}))...> powers{};
    return std::make_pair(constants, powers);
}

} // namespace detail

template <std::size_t I, class... Ps>
constexpr auto partials_with_multipliers(PowersList<Ps...>)
{
    static_assert(I < sizeof...(Ps), "Partial derivative index out of bounds");
    constexpr auto indices = typename detail::FilterZeroPowers<I, Ps...>::seq_type{};
    constexpr auto filtered = typename detail::FilterZeroPowers<I, Ps...>::plist_type{};
    constexpr auto constants_and_powers = detail::do_partials<I>(filtered);
    return std::tuple(indices, constants_and_powers.first, constants_and_powers.second);
}

} // namespace Polynomials

} // namespace Galerkin

#endif // POLYNOMIAL_POWERS_HPP

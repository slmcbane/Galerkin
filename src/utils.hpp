/*
 * Copyright (c) 2019, Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef UTILS_HPP
#define UTILS_HPP

/*!
 * @file utils.hpp
 * @brief Supporting utilities for the Galerkin library.
 * @author Sean McBane <sean.mcbane@protonmail.com>
 */

#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

/*!
 * @brief All library functionality is contained in the Galerkin namespace.
 */
namespace Galerkin
{

namespace
{
    template <class T, int N>
    struct ntuple_impl
    {
        typedef decltype(
            std::tuple_cat(
                std::make_tuple(std::declval<T>()),
                std::declval<typename ntuple_impl<T, N-1>::type>()
            )
        ) type;
    };

    template<class T>
    struct ntuple_impl<T, 0>
    {
        typedef std::tuple<> type;
    };
}

/// ntuple is a convenient alias for a tuple of N elements of the same type.
template <class T, int N>
using ntuple = typename ntuple_impl<T, N>::type;

/// Template constant yields 1 cast to the given type.
template <class T>
constexpr T one = static_cast<T>(1);

/// Template constant yields 0 cast to the given type.
template <class T>
constexpr T zero = static_cast<T>(0);

/*!
 * @brief Do an arbitrary reduction on possibly compile-time expressions.
 * 
 * This performs the reduction
 * 
 *     auto y = x;
 *     for (auto i = BEGIN; i < END; i += STEP)
 *     {
 *         y = c(y, f(i));
 *     }
 *     return y;
 * 
 * The above loop must be "type-stable", but this facility performs the operation
 * in a recursive manner so the type of `y` might change every loop iteration.
 * In addition, `f` receives the "loop index" as an `integral_constant` so that
 * the function argument can be used to specify a template parameter (e.g. `std::get`).
 * 
 * Useful for doing sums and products with compile-time types like `Rational` and
 * `Multinomial`.
 * 
 * @tparam BEGIN The initial value for the reduction loop
 * @tparam END The end value for the reduction loop; range is `[BEGIN, END)`.
 * @tparam STEP The increment between evaluations.
 * @param f A generic callable object that accepts a `std::integral_constant`.
 * @param x An initial value for the reduction - e.g. `zero<T>` for a sum or
 * `one<T>` for a product.
 * @param c A callable taking two objects `x` and `y` and returning a single
 * value; for example, `operator+` to perform a summation.
 */
template <auto BEGIN, auto END, auto STEP, class F, class COMB, class T>
constexpr auto static_reduce(F&& f, T x, COMB&& c)
{
    static_assert(BEGIN <= END);

    if constexpr (BEGIN >= END)
    {
        return x;
    }
    else
    {
        auto y = std::forward<F>(f)(std::integral_constant<decltype(BEGIN), BEGIN>());
        return static_reduce<BEGIN+STEP, END, STEP>(
            std::forward<F>(f), std::forward<COMB>(c)(x, y), std::forward<COMB>(c)
        );
    }
}

template <auto BEGIN, auto END, class F, class T, auto STEP=1>
constexpr auto static_sum(F&& f, T x)
{
    constexpr auto comb = [](auto x, auto y) { return x + y; };
    return static_reduce<BEGIN, END, STEP>(std::forward<F>(f), x, comb);
}

// This is a type-level list. It is intended to hold values that represent some
// sort of mathematical constant and have a weak ordering with <=.
template <class... Types>
struct typeconst_list;

template <class... Types>
constexpr auto make_list(Types...)
{
    return typeconst_list<Types...>();
}

template<>
struct typeconst_list<>
{
    static constexpr auto sorted()
    {
        return typeconst_list<>();
    }

    template <class... Types>
    static constexpr auto append(typeconst_list<Types...>)
    {
        return typeconst_list<Types...>();
    }

    static constexpr auto count() { return 0UL; }

    template <class F>
    static constexpr auto map([[maybe_unused]] F&& f) { return typeconst_list<>(); }

    static constexpr auto unique() { return typeconst_list<>(); }
};

template <class T, class... Types>
struct typeconst_list<T, Types...>
{
    static constexpr auto sorted()
    {
        // base case, only one type in the list.
        if constexpr (sizeof...(Types) == 0)
        {
            return typeconst_list<T>();
        }
        else
        {
            constexpr auto sorted_tail = tail().sorted();
            if constexpr(head() <= sorted_tail.head())
            {
                return typeconst_list<T>().append(sorted_tail);
            }
            else
            {
                constexpr auto first_two = make_list(sorted_tail.head(), T());
                return first_two.append(sorted_tail.tail()).sorted();
            }
        }
    }

    static constexpr auto unique()
    {
        if constexpr (sizeof...(Types) == 0)
        {
            return typeconst_list<T>();
        }
        else
        {
            if constexpr (T() == tail().head())
            {
                return tail().unique();
            }
            else
            {
                return make_list(T()).append(tail().unique());
            }
        }
    }

    template <class... OtherTypes>
    static constexpr auto append(typeconst_list<OtherTypes...>)
    {
        return typeconst_list<T, Types..., OtherTypes...>();
    }

    static constexpr auto head() { return T(); }

    static constexpr auto tail() { return typeconst_list<Types...>(); }

    static constexpr auto count() { return 1 + sizeof...(Types); }

    template <class F>
    static constexpr auto map(F&& f)
    {
        if constexpr (count() == 1)
        {
            return make_list(std::forward<F>(f)(T()));
        }
        else
        {
            return make_list(std::forward<F>(f)(T())).append(tail().map(std::forward<F>(f)));
        }
    }
};

template <auto I, class... Types>
constexpr auto get(typeconst_list<Types...> lst)
{
    static_assert(I >= 0 && I < lst.count(), "Out of bounds access to list");
    if constexpr (I == 0)
    {
        return lst.head();
    }
    else
    {
        return get<I-1>(lst.tail());
    }
}

} /* namespace Galerkin */

#endif /* UTILS_HPP */
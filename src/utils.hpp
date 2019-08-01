/*
 * Copyright (c) 2019, Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef UTILS_HPP
#define UTILS_HPP

#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Galerkin
{

using std::tuple;
using std::tuple_cat;
using std::true_type;
using std::false_type;

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

template <class T, int N>
using ntuple = typename ntuple_impl<T, N>::type;

template <class T>
constexpr T one() { return static_cast<T>(1); }

template <class T>
constexpr T zero() { return static_cast<T>(0); }

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
    static constexpr auto map(F&& f) { return typeconst_list<>(); }

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
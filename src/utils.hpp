#ifndef UTILS_HPP
#define UTILS_HPP

#include <functional>
#include <tuple>
#include <utility>

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

    if constexpr (BEGIN == END)
    {
        return x;
    }
    else
    {
        auto y = std::forward<F>(f)(std::integral_constant<decltype(BEGIN), BEGIN>());
        return static_reduce<BEGIN+1, END, STEP>(
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

} /* namespace Galerkin */

#endif /* UTILS_HPP */
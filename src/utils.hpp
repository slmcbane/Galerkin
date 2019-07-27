#ifndef UTILS_HPP
#define UTILS_HPP

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

} /* namespace Galerkin */

#endif /* UTILS_HPP */
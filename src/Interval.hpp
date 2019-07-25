#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include "utils.hpp"
#include "Interval-detail.hpp"

namespace Galerkin
{

template <class T, int ORDER>
struct IntervalElement
{
    constexpr IntervalElement() noexcept = default;

    template <int I>
    constexpr static auto basis_fn() noexcept
    {
        return Detail::interval_basis<I, T, ORDER>();
    }

    template <int I>
    constexpr static auto basis_fn_gradient(ntuple<T, ORDER+1> nodes) noexcept
    {
        return Detail::interval_basis_gradient<I, T, ORDER>();
    }

    constexpr static auto jacobian_det(ntuple<T, ORDER+1> nodes) noexcept
    {
        return Detail::interval_jacobian<T, ORDER>(nodes);
    }
};


} /* namespace Galerkin */

#endif /* INTERVAL_HPP */
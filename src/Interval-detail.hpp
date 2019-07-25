#ifndef INTERVAL_DETAIL_HPP
#define INTERVAL_DETAIL_HPP

namespace Galerkin
{

namespace Detail
{
    template <int I, class T, int ORDER>
    struct interval_basis;

    template <int I, class T>
    struct interval_basis<I, T, 1>
    {
        static_assert(I == 0 || I == 1, "Invalid node index for order = 1");

        constexpr interval_basis() noexcept = default;

        constexpr T operator()(T xi) const noexcept
        {
            const T a(-1), b(1);
            if constexpr (I == 0)
            {
                return (b - xi) / (b - a);
            }
            else
            {
                return (xi - a) / (b - a);
            }
        }
    };
    
    template <class T, int ORDER>
    struct interval_jacobian;

    template <class T>
    struct interval_jacobian<T, 1>
    {
        const ntuple<T, 2> nodes;

        constexpr interval_jacobian(ntuple<T, 2> n) noexcept : nodes(n) {}

        constexpr T operator()(T xi) const noexcept
        {
            T a, b;
            std::tie(a, b) = nodes;
            return (b - a) / 2;
        }
    };

    template <class T, int ORDER>
    struct interval_basis_gradient;

    template <class T>
    struct interval_basis_gradient<T, 1>
    {
        const ntuple<T, 2> nodes;

        constexpr interval_basis_gradient(ntuple<T, 2> n) noexcept : nodes(n) {}

        constexpr T operator()(T xi) const noexcept
        {
            T a, b;
            std::tie(a, b) = nodes;
            return 2 / (b - a);
        }
    };

} /* namespace Detail */

} /* namespace Detail */

#endif /* INTERVAL_DETAIL_HPP */
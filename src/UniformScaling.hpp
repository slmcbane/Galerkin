/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef UNIFORMSCALING_HPP
#define UNIFORMSCALING_HPP

#include "TransformBase.hpp"
#include "Quadrature.hpp"
#include "Rationals.hpp"
#include "FunctionBase.hpp"

#include <array>

namespace Galerkin
{

namespace Transforms
{

template <class T, auto N>
class UniformScaling : public TransformBase<N, UniformScaling<T, N>>
{
public:
    constexpr UniformScaling(T scale, const std::array<T, N> &translation) noexcept :
        m_scaling{scale}, m_trans(translation)
    {}

    /*!
     * @brief Call operator; A is `array-like` (support indexing with []).
     */
    template <class A>
    constexpr auto operator()(const A& xi) const noexcept
    {
        return static_reduce<0, N, 1>(
            [&](auto I)
            {
                return m_scaling * std::get<I()>(xi) + m_trans[I()];
            },
            std::tuple<>(),
            [](auto tup, auto v) { return std::tuple_cat(tup, std::tuple(v)); }
        );
    }

    /*!
     * @brief The Jacobian determinant is just scaling ^ N.
     */
    constexpr auto detJ() const noexcept
    {
        T dj = m_scaling;
        for (auto i = one<decltype(N)>; i < N; ++i)
        {
            dj *= m_scaling;
        }
        return Functions::ConstantFunction(dj);
    }

    /*!
     * @brief Integrate a function over the reference interval.
     */
    template <int I, class F>
    constexpr auto quadrature(const F &f) const noexcept
    {
        constexpr auto npoints = (I + 1) / 2 + (I - 1) % 2;
        if constexpr (npoints <= 0)
        {
            return Quadrature::integrate(f, Quadrature::legendre_rule<T, 1>);
        }
        else
        {
            return Quadrature::integrate(f, Quadrature::legendre_rule<T, npoints>);
        }
    }

    template <int I, int J>
    constexpr auto jacobian() const noexcept
    {
        static_assert(I < N && J < N, "Out of bounds index for jacobian");
        if constexpr(I == J)
        {
            return Functions::ConstantFunction(m_scaling);
        }
        else
        {
            return Functions::ConstantFunction(Rationals::rational<0>);
        }
    }

    template <int I, int J>
    constexpr auto inv_jacobian() const noexcept
    {
        static_assert(I < N && J < N, "Out of bounds index for jacobian");
        if constexpr (I == J)
        {
            return Functions::ConstantFunction(1 / m_scaling);
        }
        else
        {
            return Functions::ConstantFunction(Rationals::rational<0>);
        }
    }

private:
    T m_scaling;
    std::array<T, N> m_trans;
};

/********************************************************************************
 * Tests for uniform scaling.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

#include "Metanomials.hpp"

TEST_CASE("[Galerkin::Transforms] Test uniform scaling transformation")
{

SUBCASE("A one-dimensional uniform scaling transformation with no volume change")
{
    constexpr auto transform = UniformScaling(1.0, std::array<double, 1>{1.0});
    REQUIRE(transform.detJ()(0.1) == doctest::Approx(1.0));
    REQUIRE(transform.detJ()(-0.5) == doctest::Approx(1.0));

    REQUIRE(std::get<0>(transform(std::array<double, 1>{0.0})) == doctest::Approx(1.0));
    REQUIRE(std::get<0>(transform(std::array<double, 1>{-1.0})) == doctest::Approx(0.0));
    REQUIRE(std::get<0>(transform(std::array<int, 1>{1})) == doctest::Approx(2.0));

    // Check jacobian elements.
    REQUIRE(transform.jacobian<0, 0>()(std::array<double, 1>{0.0}) == 1.0);
    REQUIRE(transform.inv_jacobian<0, 0>()(std::array<double, 1>{0.2}) == 1.0);

    // Check partial derivatives of a function.
    REQUIRE(transform.partial<0>(Metanomials::metanomial(
        Metanomials::term(Rationals::rational<2>, Metanomials::powers(intgr_constant<2>))
        ))(std::tuple(0.2)) == doctest::Approx(0.8));

    // Check integration of a function over the interval.
    REQUIRE(transform.integrate<3>(
        Metanomials::metanomial(
            Metanomials::term(Rationals::rational<1>, Metanomials::powers(intgr_constant<3>))
        )
    ) == doctest::Approx(0.0));

    // Check integrating the partial derivative.
    REQUIRE(
        transform.integrate<2>(
            transform.partial<0>(Metanomials::metanomial(
                Metanomials::term(Rationals::rational<1>, Metanomials::powers(intgr_constant<3>))
            ))) ==
        doctest::Approx(
            transform.integrate<2>(
                Metanomials::metanomial(
                    Metanomials::term(Rationals::rational<3>, Metanomials::powers(intgr_constant<2>))
                )
            )
        )
    );
}

SUBCASE("A one dimensional uniform transformation with volume change")
{
    constexpr auto transform = UniformScaling(0.5, std::array<double, 1>{0.25});

    REQUIRE(transform.detJ()(0.1) == doctest::Approx(0.5));
    REQUIRE(transform.detJ()(-0.5) == doctest::Approx(0.5));

    REQUIRE(std::get<0>(transform(std::array<double, 1>{0.0})) == doctest::Approx(0.25));
    REQUIRE(std::get<0>(transform(std::array<double, 1>{-1.0})) == doctest::Approx(-0.25));
    // Call with another argument type.
    REQUIRE(std::get<0>(transform(std::tuple<int>(1))) == doctest::Approx(0.75));

    // Check jacobian elements.
    REQUIRE(transform.jacobian<0, 0>()(std::array<double, 1>{0.0}) == 0.5);
    REQUIRE(transform.inv_jacobian<0, 0>()(std::array<double, 1>{0.2}) == doctest::Approx(2.0));

    // Check partial derivatives of a function.
    REQUIRE(transform.partial<0>(Metanomials::metanomial(
        Metanomials::term(Rationals::rational<2>, Metanomials::powers(intgr_constant<2>))
        ))(std::tuple(0.2)) == doctest::Approx(1.6));

    // Check integration of a function over the interval.
    REQUIRE(transform.integrate<3>(
        Metanomials::metanomial(
            Metanomials::term(Rationals::rational<1>, Metanomials::powers(intgr_constant<3>)),
            Metanomials::term(Rationals::rational<-1>, Metanomials::powers(intgr_constant<2>))
        )
    ) == doctest::Approx(-1.0 / 3));

    // Check integrating the partial derivative.
    REQUIRE(
        transform.integrate<2>(
            transform.partial<0>(Metanomials::metanomial(
                Metanomials::term(Rationals::rational<1>, Metanomials::powers(intgr_constant<3>))
            ))) ==
        doctest::Approx(
            transform.integrate<2>(
                Metanomials::metanomial(
                    Metanomials::term(Rationals::rational<6>, Metanomials::powers(intgr_constant<2>))
                )
            )
        )
    );

    REQUIRE(
        transform.integrate<2>(
            transform.partial<0>(Metanomials::metanomial(
                Metanomials::term(Rationals::rational<1>, Metanomials::powers(intgr_constant<3>))
            ))) == doctest::Approx(2.0)
    );

    REQUIRE(
        transform.integrate<2>(
            Metanomials::metanomial(
                    Metanomials::term(Rationals::rational<3>, Metanomials::powers(intgr_constant<2>))
                )
        ) == doctest::Approx(1.0)
    );
}

} /* TEST_CASE */

#endif

} // namespace Transforms

} // namespace Galerkin

#endif /* UNIFORMSCALING_HPP */

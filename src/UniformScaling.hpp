/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef UNIFORMSCALING_HPP
#define UNIFORMSCALING_HPP

#include "TransformBase.hpp"

#include <array>

namespace Galerkin
{

namespace Transforms
{

template <class T, auto N>
class UniformScaling : public TransformBase<N, UniformScaling<T, N>>
{
    /*!
     * @brief Call operator; A is `array-like` (support indexing with []).
     */
    template <class A>
    constexpr auto operator()(const A& xi) const noexcept
    {
        std::array<T, N> x{};
        for (auto i = zero<decltype(N)>; i < N; ++i)
        {
            x[i] = xi[i] * m_scaling + m_trans[i];
        }
        return x;
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
    auto transform = uniform_transformation(1.0, std::array<double, 1>{1.0});
    REQUIRE(transform.detJ()(0.1) == doctest::Approx(1.0));
    REQUIRE(transform.detJ()(-0.5) == doctest::Approx(1.0));

    REQUIRE(transform(std::array<double, 1>{0.0})[0] == doctest::Approx(1.0));
    REQUIRE(transform(std::array<double, 1>{-1.0})[0] == doctest::Approx(0.0));
    REQUIRE(transform(std::array<int, 1>{1})[0] == doctest::Approx(2.0));

    // Check jacobian elements.
    REQUIRE(transform.jacobian<0, 0>(std::array<double, 1>{0.0}) == 1.0);
    REQUIRE(transform.inv_jacobian<0, 0>(std::array<double, 1>{0.2}) == 1.0);

    // Check partial derivatives of a function.
    REQUIRE(transform.partial<0>(Metanomials::metanomial(
        Metanomials::term(Rationals::rational<2>, Metanomials::powers(intgr_constant<2>))
        ))(std::array<double, 1>{0.2}) == doctest::Approx(0.8));

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

} /* TEST_CASE */

#endif

} // namespace Transforms

} // namespace Galerkin

#endif /* UNIFORMSCALING_HPP */

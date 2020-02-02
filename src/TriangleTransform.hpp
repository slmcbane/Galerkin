/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef TRIANGLETRANSFORM_HPP
#define TRIANGLETRANSFORM_HPP

/*!
 * @file TriangleTransform.hpp
 * @brief Implementation of transformation from reference to instantiated triangle domain.
 */

#include "Rationals.hpp"
#include "TransformBase.hpp"

#include <array>
#include <type_traits>

namespace Galerkin
{

namespace Transforms
{

template <class T>
class TriangleTransform : public TransformBase<2, TriangleTransform<T>>
{
public:
    template <class P1, class P2, class P3>
    constexpr TriangleTransform(const P1 &p1, const P2 &p2, const P3 &p3) noexcept :
        m_coeffs{ 0 }
    {
        auto x1 = get<0>(p1); auto y1 = get<1>(p1);
        auto x2 = get<0>(p2); auto y2 = get<1>(p2);
        auto x3 = get<0>(p3); auto y3 = get<1>(p3);

        m_coeffs[0] = x3 - x1;
        m_coeffs[1] = x2 - x1;
        m_coeffs[2] = x2 + x3;
        m_coeffs[3] = y3 - y1;
        m_coeffs[4] = y2 - y1;
        m_coeffs[5] = y2 + y3;
        for (auto &coeff: m_coeffs)
        {
            coeff /= 2;
        }
    }

    template <class Arg>
    constexpr auto operator()(const Arg &arg) const noexcept
    {
        const auto xi = get<0>(arg);
        const auto eta = get<1>(arg);
        const auto x = m_coeffs[0] * xi + m_coeffs[1] * eta + m_coeffs[2];
        const auto y = m_coeffs[3] * xi + m_coeffs[4] * eta + m_coeffs[5];
        return std::array<std::remove_cv_t<decltype(x)>, 2>{ x, y };
    }

private:
    std::array<T, 6> m_coeffs;
};

template <class... Ps>
constexpr auto triangle_transform(const Ps &...ps) noexcept
{
    static_assert(sizeof...(Ps) == 3);
    return TriangleTransform<double>(ps...);
}

/********************************************************************************
 * Test block for triangle transform
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Transforms] Test mapping of a triangle")
{

constexpr auto transform = triangle_transform(
    std::tuple(0.5, 0.5), std::tuple(0.0, 1.0), std::tuple(1.0, 1.0)
);

SUBCASE("Check transformations of points")
{
    auto pt = transform(std::tuple(-1.0, -1.0));
    REQUIRE(get<0>(pt) == doctest::Approx(0.5));
    REQUIRE(get<1>(pt) == doctest::Approx(0.5));

    pt = transform(std::tuple(-0.5, -0.5));
    REQUIRE(get<0>(pt) == doctest::Approx(0.5));
    REQUIRE(get<1>(pt) == doctest::Approx(0.75));

    pt = transform(std::tuple(-0.5, 0.5));
    REQUIRE(get<0>(pt) == doctest::Approx(0.25));
    REQUIRE(get<1>(pt) == doctest::Approx(1.0));
} // SUBCASE

} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace Transforms

} // namespace Galerkin

#endif // TRIANGLETRANSFORM_HPP

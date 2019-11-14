#ifndef TRANSFORMS_HPP
#define TRANSFORMS_HPP

#include <Eigen/Core>

namespace Galerkin
{

/*!
 * @brief Namespace for functionality related to coordinate transformations.
 */
namespace Transforms
{

template <auto N, class Derived>
class TransformBase
{
public:

    template <class T>
    auto operator()(T... xs) const noexcept
    {
        static_assert(sizeof...(xs) == N, "Wrong number of arguments for transform dimensionality");
        return static_cast<const Derived*>(this)->apply(Eigen::Matrix<T, N, 1>(xs...));
    }

    template <class T>
    auto operator()(const Eigen::Matrix<T, N, 1>& x) const noexcept
    {
        return static_cast<const Derived*>(this)->apply(x);
    }

    
private:
    TransformBase() = default;
    friend Derived;
};

/********************************************************************************
 * Tests for uniform and diagonal scalings.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

#include "Multinomials.hpp"

TEST_CASE("[Galerkin::Transforms] Test uniform transformation")
{

SUBCASE("A one-dimensional uniform transformation with no volume change")
{
    auto transform = uniform_transformation(std::tuple(0.0, 2.0));
    REQUIRE(transform.detJ(0.1) == doctest::Approx(1.0));
    REQUIRE(transform.detJ(1.1) == doctest::Approx(1.0));
    REQUIRE(transform(0.0)(0) == doctest::Approx(1.0));
    // Check the vector form.
    REQUIRE(transform(0.0).isApprox(Eigen::Matrix<double, 1, 1>(1.0)));

    REQUIRE(transform(-1.0)(0) == doctest::Approx(0.0));

    REQUIRE(transform(1.0) == doctest::Approx(2.0));
    REQUIRE(transform(1.0).isApprox(Eigen::Matrix<double, 1, 1>(2.0)));
    // Should also accept an Eigen::Vector with the appropriate size.
    REQUIRE(transform(Eigen::Matrix<double, 1, 1>(0.0)).isApprox(Eigen::Matrix<double, 1, 1>(1.0)));

    // Check jacobian elements.
    REQUIRE(transform.jacobian<0, 0>(0.0) == 1.0);
    REQUIRE(transform.inv_jacobian<0, 0>(0.2) == 1.0);
    REQUIRE(transform.inv_jacobian<0, 0>(Eigen::Matrix<double, 1, 1>(1.3)) == 1.0);

    // Check partial derivatives of a function.
    REQUIRE(transform.partial<0>(Multinomials::multinomial(
        Multinomials::term(Rationals::rational<2>, Multinomials::powers(intgr_constant<2>))
        ))(0.2) == doctest::Approx(0.8));
    REQUIRE(transform.partial<0>(Multinomials::multinomial(
        Multinomials::term(Rationals::rational<2>, Multinomials::powers(intgr_constant<2>))
        ))(Eigen::Matrix<double, 1, 1>(0.2)) == doctest::Approx(0.8));

    // Check integration of a function over the interval.
    REQUIRE(transform.integrate<3>(
        Multinomials::multinomial(
            Multinomials::term(Rationals::rational<1>, Multinomials::powers(intgr_constant<3>))
        )
    ) == doctest::Approx(0.0));

    // Check integrating the partial derivative.
    REQUIRE(
        transform.integrate<2>(
            transform.partial<0>(Multinomials::multinomial(
                Multinomials::term(Rationals::rational<1>, Multinomials::powers(intgr_constant<3>))
            ))) ==
        doctest::Approx(
            transform.integrate<2>(
                Multinomials::multinomial(
                    Multinomials::term(Rationals::rational<3>, Multinomials::powers(intgr_constant<2>))
                )
            )
        )
    );
}

} /* TEST_CASE */

#endif

} /* namespace Transforms */

} /* namespace Galerkin */

#endif /* TRANSFORMS_HPP */

/*
 * Copyright (c) 2019, The University of Texas at Austin & Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef ELEMENTS_HPP
#define ELEMENTS_HPP

#include "MetaLinAlg.hpp"
#include "Multinomials.hpp"
#include "utils.hpp"

#include <tuple>

namespace Galerkin
{

/// Functionality related to individual finite elements.
namespace Elements
{

/*!
 * @brief Represents the form of a multinomial shape function.
 * 
 * The `ShapeFunctionForm` is just an empty variadic template struct, where the
 * template parameters are types of `Powers` classes representing the terms of a
 * multinomial. For example, the following:
 * 
 *     constexpr ShapeFunctionForm<Powers<0, 0>, Powers<0, 1>, Powers<1, 0>, Powers<1, 1>> form{};
 * 
 * represents the form of a bilinear shape function `phi(x, y) = a*x*y + b*x + c*y + d`.
 */
template <class... Powers>
struct ShapeFunctionForm : public typeconst_list<Powers...>
{
};

// Helper to make sure arguments to make_form are all Powers objects.
namespace
{

template <class P, class... Ps>
constexpr bool check_powers() noexcept
{
    if constexpr (sizeof...(Ps) == 0)
    {
        return Multinomials::is_powers<P>;
    }
    else
    {
        return Multinomials::is_powers<P> && check_powers<Ps...>();
    }
}

template <class... Ps>
constexpr auto to_form(typeconst_list<Ps...>) noexcept
{
    return ShapeFunctionForm<Ps...>();
}

} // namespace

/*!
 * @brief Construct a `ShapeFunctionForm` from variadic list of `Powers` objects
 * 
 * The initialized form has duplicate `Powers` objects merged into one and terms
 * sorted in ascending order. Arguments are checked to make sure that they are
 * actually `Powers`.
 */
template <class... Powers>
constexpr auto make_form(Powers...) noexcept
{
    static_assert(check_powers<Powers...>(), "All arguments to make_form should be 'Powers' objects");
    return to_form(ShapeFunctionForm<Powers...>().sorted().unique());
}

/*!
 * @brief Struct representing a "control point", or vertex of an element.
 * 
 * A `ControlPoint` is one of the vertices of an abstract reference element. I.E.,
 * for a classical quadrilateral element, the control points are (probably)
 * `(-1, -1)`, `(-1, 1)`, `(1, 1)`, `(1, -1)`. This struct is only a lightweight
 * tag object but it does have the additional functionality of conversion to a
 * tuple, for use in evaluating polynomial terms at the point.
 */
template <class... Coords>
struct ControlPoint
{
    /// Convert `ControlPoint` to tuple for use as a function argument.
    static constexpr auto to_tuple() noexcept
    {
        return std::tuple(Coords()...);
    }
};

// Helper function to check types of coordinate variables.
namespace
{
template <class Coord, class... Coords>
constexpr bool check_coords() noexcept
{
    constexpr bool is_coord = Rationals::is_rational<Coord> || is_intgr_constant<Coord>;
    if constexpr (sizeof...(Coords) == 0)
    {
        return is_coord;
    }
    else
    {
        return is_coord && check_coords<Coords...>();
    }
}

} // namespace

/*!
 * @brief Construct a control point from list of coordinates.
 * 
 * Each coordinate argument to this function should be either a `Rational` or an
 * `integral_constant`. This is checked. The returned value is an N-dimensional
 * control point, where `N` is the number of coordinates given.
 * 
 * @param[in] coords... List of coordinates of the point.
 */
template <class... Coords>
constexpr auto control_point(Coords...)
{
    static_assert(check_coords<Coords...>(), "All coordinates should be rationals or integral_constants");
    return ControlPoint<Coords...>();
}

// Here ends boilerplate for the DSL and begins the actual implementation of
// deriving shape functions.
namespace
{

template <class... Powers, class Pt>
constexpr auto compute_power_values(typeconst_list<Powers...>, Pt) noexcept
{
    return make_list(Multinomials::raise(Pt::to_tuple(), Powers())...);
}

template <class... Powers, class... Pts>
constexpr auto build_terms_matrix(ShapeFunctionForm<Powers...>, typeconst_list<Pts...>) noexcept
{
    return static_reduce<0, sizeof...(Pts), 1>(
        [](auto I) {
            return compute_power_values(typeconst_list<Powers...>(), get<I()>(typeconst_list<Pts...>()));
        },
        make_list(),
        [](auto L, auto x) { return L.append(make_list(x)); });
}

template <class... Coeffs, class... Powers>
constexpr auto multiply_coeffs(typeconst_list<Coeffs...>, ShapeFunctionForm<Powers...>) noexcept
{
    return Multinomials::multinomial(Multinomials::term(Coeffs(), Powers())...);
}

} // namespace

/*!
 * @brief Given nodes (control points) of an element, derive shape functions.
 * 
 * This function derives multinomial shape functions on an element given the
 * locations of its nodes. The functions are derived by setting them equal to 1
 * at nodes, and equal to 0 on the others. The length of `Powers...` must equal
 * the length of `Pts...` to obtain a well-posed system.
 */
template <class... Powers, class... Pts>
constexpr auto derive_shape_functions(ShapeFunctionForm<Powers...>, typeconst_list<Pts...>) noexcept
{
    static_assert(sizeof...(Powers) == sizeof...(Pts), "Ill-posed system to derive shape functions");
    constexpr auto terms_matrix = build_terms_matrix(
        ShapeFunctionForm<Powers...>(), typeconst_list<Pts...>());

    return static_reduce<0, sizeof...(Pts), 1>(
        [=](auto I) {
            constexpr auto coeffs = MetaLinAlg::linear_solve(terms_matrix,
                                                             MetaLinAlg::canonical<I(), sizeof...(Pts)>());
            return multiply_coeffs(coeffs, ShapeFunctionForm<Powers...>());
        },
        typeconst_list<>(),
        [](auto L, auto x) { return L.append(make_list(x)); });
}

/********************************************************************************
 * Test derivation of shape functions given "control points" and a form for the
 * shape function.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using namespace Multinomials;
using namespace Rationals;

TEST_CASE("[Galerkin::Elements] Deriving multinomial shape functions")
{
    SUBCASE("Test for a one-dimensional element")
    {
        // Shape function has the form ax + b.
        constexpr auto form = make_form(
            powers(intgr_constant<1>),
            powers(intgr_constant<0>));

        // Control points are -1, 1.
        constexpr auto control_points = make_list(
            control_point(rational<-1>),
            control_point(rational<1>));

        constexpr auto fns = derive_shape_functions(form, control_points);

        REQUIRE(get<0>(fns) ==
                multinomial(
                    term(-rational<1, 2>, powers(intgr_constant<1>)),
                    term(rational<1, 2>, powers(intgr_constant<0>))));
        REQUIRE(get<1>(fns) ==
                multinomial(
                    term(rational<1, 2>, powers(intgr_constant<1>)),
                    term(rational<1, 2>, powers(intgr_constant<0>))));
    }

    SUBCASE("Test for a quadratic one-dimensional element")
    {
        constexpr auto form = make_form(
            powers(intgr_constant<2>),
            powers(intgr_constant<1>),
            powers(intgr_constant<0>));

        constexpr auto control_points = make_list(
            control_point(rational<-1>),
            control_point(rational<0>),
            control_point(rational<1>));
        
        constexpr auto fns = derive_shape_functions(form, control_points);

        REQUIRE(get<0>(fns) ==
            multinomial(
                term(rational<1, 2>, powers(intgr_constant<2>)),
                term(-rational<1, 2>, powers(intgr_constant<1>))
            ));

        REQUIRE(get<1>(fns) ==
            multinomial(
                term(-rational<1>, powers(intgr_constant<2>)),
                term(rational<1>, powers(intgr_constant<0>))
            ));

        REQUIRE(get<2>(fns) ==
            multinomial(
                term(rational<1, 2>, powers(intgr_constant<2>)),
                term(rational<1, 2>, powers(intgr_constant<1>))
            ));
    }
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 * End shape function API test.
 *******************************************************************************/
} // namespace Elements

} /* namespace Galerkin */

#endif /* ELEMENTS_HPP */
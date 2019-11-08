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

// Here ends boilerplate for the DSL and begins the actual implementation of
// deriving shape functions.
namespace
{

template <class... Powers, class... Constraints>
constexpr auto build_terms_matrix(ShapeFunctionForm<Powers...>, typeconst_list<Constraints...>) noexcept
{
    return typeconst_list<Constraints...>().map(
        [](auto constraint)
        {
            return make_list(constraint(Powers())...);
        }
    );
}

template <class... Coeffs, class... Powers>
constexpr auto multiply_coeffs(typeconst_list<Coeffs...>, ShapeFunctionForm<Powers...>) noexcept
{
    return Multinomials::multinomial(Multinomials::term(Coeffs(), Powers())...);
}

} // namespace

/*!
 * @brief Given degrees of freedom for an element, derive shape functions.
 * 
 * This function derives multinomial shape functions on an element given the
 * degrees of freedom on the element in a functional form. The first argument is
 * a `ShapeFunctionForm` specifying the powers of the multinomial form, and the
 * second is a `typeconst_list` of constraints. These take the form of a function
 * object accepting a `Powers` object (see `Multinomials.hpp` for interface) and
 * returning a number as a `Rational`. Functors for the most common cases are
 * provided - see `evaluate_at` and `partial_at`.
 */
template <class... Powers, class... Constraints>
constexpr auto derive_shape_functions(ShapeFunctionForm<Powers...>, typeconst_list<Constraints...>) noexcept
{
    static_assert(sizeof...(Powers) == sizeof...(Constraints), "Ill-posed system to derive shape functions");
    constexpr auto terms_matrix = build_terms_matrix(
        ShapeFunctionForm<Powers...>(), typeconst_list<Constraints...>());

    return static_reduce<0, sizeof...(Constraints), 1>(
        [=](auto I) {
            constexpr auto coeffs = MetaLinAlg::linear_solve(terms_matrix,
                                                             MetaLinAlg::canonical<I(), sizeof...(Constraints)>());
            return multiply_coeffs(coeffs, ShapeFunctionForm<Powers...>());
        },
        typeconst_list<>(),
        [](auto L, auto x) { return L.append(make_list(x)); });
}

template <class... Ns>
struct EvaluateAt
{
    template <class Powers>
    constexpr auto operator()(Powers) const noexcept
    {
        return Multinomials::raise(std::tuple(Ns()...), Powers());
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

template <class... Ns>
constexpr auto evaluate_at(Ns...) noexcept
{
    static_assert(check_coords<Ns...>(), "All coordinates should be rationals or integral_constants");
    return EvaluateAt<Ns...>();
}

template <class CoordList, auto I, auto... Is>
class PartialAt
{
public:
    template <class Powers>
    constexpr auto operator()(Powers) const noexcept
    {
        constexpr auto t = take_partials<I, Is...>(Powers());
        return t(instantiate_tuple(CoordList()));
    }

private:

    template <auto J, auto... Js, class Powers>
    static constexpr auto take_partials(Powers) noexcept
    {
        constexpr auto first_partial = Multinomials::partial<J>(
            Multinomials::term(Rationals::rational<1>, Powers())
        );

        if constexpr (sizeof...(Js) == 0)
        {
            return first_partial;
        }
        else
        {
            return first_partial.coeff() * take_partials<Js...>(first_partial.powers());
        }
    }

    template <class... Coords>
    static constexpr auto instantiate_tuple(typeconst_list<Coords...>) noexcept
    {
        return std::tuple(Coords()...);
    }
};

template <auto I, auto... Is, class... Coords>
constexpr auto partial_at(Coords...) noexcept
{
    static_assert(check_coords<Coords...>(), "All coordinates should be rationals or integral_constants");
    return PartialAt<typeconst_list<Coords...>, I, Is...>();
}

/********************************************************************************
 * Test derivation of shape functions given "control points" and a form for the
 * shape function.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using namespace Multinomials;
using namespace Rationals;

TEST_CASE("[Galerkin::Elements] Deriving one-dimensional shape functions")
{
    SUBCASE("Test for a first order element")
    {
        // Shape function has the form ax + b.
        constexpr auto form = make_form(
            powers(intgr_constant<1>),
            powers(intgr_constant<0>));

        // Constraints are the function values at -1, 1.
        constexpr auto constraints = make_list(
            evaluate_at(rational<-1>),
            evaluate_at(rational<1>));

        constexpr auto fns = derive_shape_functions(form, constraints);

        REQUIRE(get<0>(fns) ==
                multinomial(
                    term(-rational<1, 2>, powers(intgr_constant<1>)),
                    term(rational<1, 2>, powers(intgr_constant<0>))));
        REQUIRE(get<1>(fns) ==
                multinomial(
                    term(rational<1, 2>, powers(intgr_constant<1>)),
                    term(rational<1, 2>, powers(intgr_constant<0>))));
    }

    SUBCASE("Test for a second order element")
    {
        constexpr auto form = make_form(
            powers(intgr_constant<2>),
            powers(intgr_constant<1>),
            powers(intgr_constant<0>));

        constexpr auto constraints = make_list(
            evaluate_at(rational<-1>),
            evaluate_at(rational<0>),
            evaluate_at(rational<1>));
        
        constexpr auto fns = derive_shape_functions(form, constraints);

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

    SUBCASE("Test for a 3rd-order element with derivative DOF's")
    {
        constexpr auto form = make_form(
            powers(intgr_constant<3>),
            powers(intgr_constant<2>),
            powers(intgr_constant<1>),
            powers(intgr_constant<0>));

        constexpr auto constraints = make_list(
            evaluate_at(rational<-1>),
            evaluate_at(rational<1>),
            partial_at<0>(rational<-1>),
            partial_at<0>(rational<1>));

        constexpr auto fns = derive_shape_functions(form, constraints);

        REQUIRE(get<0>(fns) ==
            multinomial(
                term(rational<1, 4>, powers(intgr_constant<3>)),
                term(-rational<3, 4>, powers(intgr_constant<1>)),
                term(rational<1, 2>, powers(intgr_constant<0>))
            ));

        REQUIRE(get<1>(fns) ==
            multinomial(
                term(-rational<1, 4>, powers(intgr_constant<3>)),
                term(rational<3, 4>, powers(intgr_constant<1>)),
                term(rational<1, 2>, powers(intgr_constant<0>))
            ));

        REQUIRE(get<2>(fns) ==
            multinomial(
                term(rational<1, 4>, powers(intgr_constant<3>)),
                term(-rational<1, 4>, powers(intgr_constant<2>)),
                term(-rational<1, 4>, powers(intgr_constant<1>)),
                term(rational<1, 4>, powers(intgr_constant<0>))
            ));

        REQUIRE(get<3>(fns) ==
            multinomial(
                term(rational<1, 4>, powers(intgr_constant<3>)),
                term(rational<1, 4>, powers(intgr_constant<2>)),
                term(-rational<1, 4>, powers(intgr_constant<1>)),
                term(-rational<1, 4>, powers(intgr_constant<0>))
            ));
    }
}

TEST_CASE("[Galerkin::Elements] Deriving two-dimensional shape functions")
{
    SUBCASE("Test derivation of bilinear shape functions on a quadrilateral")
    {
        constexpr auto form = make_form(
            powers(intgr_constant<1>, intgr_constant<1>),
            powers(intgr_constant<1>, intgr_constant<0>),
            powers(intgr_constant<0>, intgr_constant<1>),
            powers(intgr_constant<0>, intgr_constant<0>)
        );

        constexpr auto constraints = make_list(
            evaluate_at(rational<-1>, rational<-1>),
            evaluate_at(rational<-1>, rational<1>),
            evaluate_at(rational<1>, rational<1>),
            evaluate_at(rational<1>, rational<-1>)
        );

        constexpr auto fns = derive_shape_functions(form, constraints);

        REQUIRE(get<0>(fns) ==
            multinomial(
                term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
                term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
                term(-rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
                term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
            ));

        REQUIRE(get<1>(fns) ==
            multinomial(
                term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
                term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
                term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
                term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
            ));

        REQUIRE(get<2>(fns) ==
            multinomial(
                term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
                term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
                term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
                term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
            ));

        REQUIRE(get<3>(fns) ==
            multinomial(
                term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
                term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
                term(-rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
                term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
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
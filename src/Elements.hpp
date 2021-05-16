/*
 * Copyright (c) 2019, The University of Texas at Austin & Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef ELEMENTS_HPP
#define ELEMENTS_HPP

/*!
 * @file Elements.hpp
 * @brief Basic functionality related to defining element types.
 */

#include "LinAlg.hpp"
#include "Polynomial.hpp"
#include "utils.hpp"

#include <tuple>

namespace Galerkin
{

namespace Elements
{

/*!
 * @brief Represents the form of a polynomial shape function.
 *
 * The `ShapeFunctionForm` is just an empty variadic template struct, where the
 * template parameters are types of `Powers` classes representing the terms of a
 * metanomial. For example, the following:
 *
 *     constexpr ShapeFunctionForm<Powers<0, 0>, Powers<0, 1>, Powers<1, 0>, Powers<1, 1>> form{};
 *
 * represents the form of a bilinear shape function `phi(x, y) = a*x*y + b*x + c*y + d`.
 */
template <class... Powers>
struct ShapeFunctionForm : public typeconst_list<Powers...>
{
    constexpr operator Polynomials::PowersList<Powers...>() const noexcept
    {
        return Polynomials::PowersList<Powers...>{};
    }
};

// Helper to make sure arguments to make_form are all Powers objects.
namespace
{

template <class... Ps>
constexpr bool check_powers() noexcept
{
    return Polynomials::detail::all_are_powers<Ps...>::value;
}

template <class... Ps>
constexpr auto powers_count(ShapeFunctionForm<Ps...>) noexcept
{
    return sizeof...(Ps);
}

} // namespace

/*!
 * @brief Convert a `typeconst_list` to a `ShapeFunctionForm`
 *
 * This function checks that all of the types `Ps...` are instantiations of
 * `Metanomials::Powers`, then returns a `ShapeFunctionForm` with the same
 * powers.
 */
template <class... Ps>
constexpr auto to_powerslist(typeconst_list<Ps...>) noexcept
{
    static_assert(check_powers<Ps...>(), "All arguments to to_form should be 'Powers' objects");
    return Polynomials::PowersList<Ps...>();
}

/*!
 * @brief Construct a `ShapeFunctionForm` from variadic list of `Metanomials::Powers` objects
 *
 * The initialized form has duplicate `Powers` objects merged into one and terms
 * sorted in ascending order. Arguments are checked to make sure that they are
 * actually `Powers`.
 */
template <class... Powers>
constexpr auto make_form(Powers...) noexcept
{
    static_assert(check_powers<Powers...>(), "All arguments to make_form should be 'Powers' objects");
    return to_powerslist(ShapeFunctionForm<Powers...>().sorted().unique());
}

/*!
 * @brief Combine two `Powers` objects
 *
 * This is the analogue of `std::tuple_cat` or `typeconst_list::append`.
 * Combine the powers in two `Powers` objects into one `Powers` objects;
 * e.g. `concatenate_powers(Powers<1>{}, Powers<1>{}) == Powers<1, 1>{}`.
 */
template <auto... Is, auto... Js>
constexpr auto concatenate_powers(Polynomials::Powers<Is...>, Polynomials::Powers<Js...>) noexcept
{
    return Polynomials::Powers<Is..., Js...>{};
}

/*!
 * @brief Returns a `ShapeFunctionForm` with combination of powers up to given maxima.
 *
 * This function is a utility to construct the common form for shape functions
 * that consists of all terms that are at most order N in a given term.
 * For example, `powers_up_to(intgr_constant<1>, intgr_constant<1>)` gives the
 * form of a bilinear shape function `f(x, y) = axy + bx + cy + d`.
 */
template <auto I, auto... Is>
constexpr auto
powers_up_to() noexcept
{
    static_assert(I >= 0);
    if constexpr (sizeof...(Is) == 0)
    {
        constexpr auto lst = static_reduce<0, I + 1, 1>(
            [](auto i) { return i; }, typeconst_list<>{},
            [](auto l, auto i) { return l.append(typeconst_list<Polynomials::Powers<i()>>{}); });
        return lst;
    }
    else
    {
        constexpr auto lst = powers_up_to<I>();
        constexpr auto tails = powers_up_to<Is...>();
        constexpr auto power_list = static_reduce<0, lst.count, 1>(
            [=](auto i) {
                constexpr auto index1 = i();
                return static_reduce<0, tails.count, 1>(
                    [=](auto j) { return concatenate_powers(get<index1>(lst), get<j()>(tails)); },
                    typeconst_list<>{},
                    [](auto l1, auto pow) { return l1.append(typeconst_list<decltype(pow)>{}); });
            },
            typeconst_list<>{}, [](auto l1, auto l2) { return l1.append(l2); });
        return power_list;
    }
}

// Here ends boilerplate for the DSL and begins the actual implementation of
// deriving shape functions.
namespace
{

template <class... Powers, class... Constraints>
constexpr auto build_terms_matrix(
    Polynomials::PowersList<Powers...>, const std::tuple<Constraints...> &constraints) noexcept
{
    LinAlg::Matrix<sizeof...(Constraints)> A;
    static_for<0, sizeof...(Constraints), 1>([&](auto I) {
        static_for<0, sizeof...(Constraints), 1>([&](auto J) {
            A = A.set_entry(
                I(), J(),
                std::get<I()>(constraints)(
                    typename std::tuple_element<J(), std::tuple<Powers...>>::type{}));
        });
    });
    return A;
}

template <class... Coeffs, class... Powers>
constexpr auto multiply_coeffs(typeconst_list<Coeffs...>, ShapeFunctionForm<Powers...>) noexcept
{
    return Polynomials::make_poly(std::tuple(Coeffs()...), Polynomials::PowersList<Powers...>{});
    // return Metanomials::metanomial(Metanomials::term(Coeffs(), Powers())...);
}

} // namespace

/*!
 * @brief Given degrees of freedom for an element, derive shape functions.
 *
 * This function derives polynomial shape functions on an element given the
 * degrees of freedom on the element in a functional form. The first argument is
 * a `ShapeFunctionForm` specifying the powers of the metanomial form, and the
 * second is a `typeconst_list` of constraints. These take the form of a function
 * object accepting a `Powers` object (see `Metanomials.hpp` for interface) and
 * returning a number as a `Rationals::Rational`. Functors for the most common cases are
 * provided - see `evaluate_at` and `partial_at`.
 *
 * @see evaluate_at
 * @see partial_at
 */
template <class... Powers, class... Constraints>
constexpr auto derive_shape_functions(
    Polynomials::PowersList<Powers...>, const std::tuple<Constraints...> &constraints) noexcept
{
    static_assert(
        sizeof...(Powers) == sizeof...(Constraints), "Ill-posed system to derive shape functions");
    const auto terms_matrix = build_terms_matrix(Polynomials::PowersList<Powers...>{}, constraints);

    using result_type = decltype(Polynomials::make_poly(
        std::declval<std::array<Rationals::Rational, sizeof...(Powers)>>(),
        Polynomials::PowersList<Powers...>{}));
    std::array<result_type, sizeof...(Constraints)> shape_functions;

    for (std::size_t i = 0; i < sizeof...(Constraints); ++i)
    {
        shape_functions[i] = Polynomials::make_poly(
            LinAlg::linear_solve(terms_matrix, LinAlg::canonical<sizeof...(Constraints)>(i)).data,
            Polynomials::PowersList<Powers...>{});
    }

    return shape_functions;
}

template <class... Powers, class... Constraints>
constexpr auto derive_shape_functions(
    typeconst_list<Powers...>, const std::tuple<Constraints...> &constraints) noexcept
{
    return derive_shape_functions(Polynomials::PowersList<Powers...>{}, constraints);
}

/*!
 * @brief A constraint functor for use in `derive_shape_functions`.
 *
 * An instance of this class, when a member of the list of constraints given to
 * `derive_shape_functions`, indicates that one of the degrees of freedom
 * constraining the system is the value of a function at the point
 * `(Ns...)` in R^n (where `n == sizeof...(Ns)`). There are no data members;
 * the point at which to evaluate must be specified using either
 * a `Rationals::Rational` or `std::integral_constant`. Construct an instance
 * of this functor using the helper function `evaluate_at`.
 *
 * @see derive_shape_functions
 * @see evaluate_at
 */
template <class... Ns>
struct EvaluateAt
{
    std::tuple<Ns...> evaluation_point;
    template <class Powers>
    constexpr auto operator()(Powers) const noexcept
    {
        return Polynomials::make_poly(
            std::tuple(Rationals::rational<1>), Polynomials::PowersList<Powers>{})(evaluation_point);
    }

    constexpr EvaluateAt(Ns &&...xs) : evaluation_point(std::forward<Ns>(xs)...) {}

    constexpr EvaluateAt(const Ns &...xs) : evaluation_point(xs...) {}
};

/*!
 * @brief Functor representing a partial derivative value as DOF
 *
 * This class is intended for use as a constraint in `derive_shape_functions`.
 * The indices `(I, Is...)` represent the variable which should be
 * differentiated with respect to, and the point at which to take the partial
 * derivative is a `typeconst_list` encoded in the type `CoordList`.
 *
 * For example, a degree of freedom which is the value of the cross derivative
 * `\frac{\partial^2 f}{\partial x \partial y}` at the origin would be encoded
 * as the type `PartialAt<typeconst_list<Rationals::Rational<0, 1>, Rationals::Rational<0, 1>>, 0, 1>`.
 *
 * In practice, construct an instance using `partial_at`.
 *
 * @see partial_at
 * @see derive_shape_functions
 */
template <std::size_t I, std::size_t... Is>
struct Partial
{
    template <class... Ts>
    struct At
    {
        constexpr At(Ts &&...xs) : evaluation_point(std::forward<Ts>(xs)...) {}

        constexpr At(const Ts &...xs) : evaluation_point(xs...) {}

        template <class Powers>
        constexpr auto operator()(Powers) const noexcept
        {
            return take_partials<I, Is...>(Polynomials::make_poly(
                std::tuple(Rationals::rational<1>), Polynomials::PowersList<Powers>{}))(evaluation_point);
        }

      private:
        template <auto J, auto... Js, class Powers>
        constexpr static auto
        take_partials(const Polynomials::Polynomial<Rationals::Rational, Powers> &p) noexcept
        {
            auto first_partial = p.template partial<J>();

            if constexpr (sizeof...(Js) == 0)
            {
                return first_partial;
            }
            else
            {
                return take_partials<Js...>(first_partial);
            }
        }
        std::tuple<Ts...> evaluation_point;
    };
};

/********************************************************************************
 * Test derivation of shape functions given "control points" and a form for the
 * shape function.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using namespace Rationals;
using namespace Polynomials;

TEST_CASE("[Elements] Test powers_up_to")
{
    constexpr auto powers = powers_up_to<1, 1>();
    REQUIRE(std::is_same_v<
            decltype(powers),
            const typeconst_list<Powers<0, 0>, Powers<0, 1>, Powers<1, 0>, Powers<1, 1>>>);
}

TEST_CASE("[Elements] Deriving shape functions")
{
    SUBCASE("Test for a first order element")
    {
        // Shape function has the form ax + b.
        constexpr auto form = PowersList<Powers<1>, Powers<0>>{};

        // Constraints are the function values at -1, 1.
        constexpr auto constraints = std::tuple(EvaluateAt(rational<-1>), EvaluateAt(rational<1>));

        constexpr auto fns = derive_shape_functions(form, constraints);

        REQUIRE(
            get<0>(fns) ==
            make_poly(std::tuple(-rational<1, 2>, rational<1, 2>), PowersList<Powers<1>, Powers<0>>{}));

        REQUIRE(
            get<1>(fns) ==
            make_poly(std::tuple(rational<1, 2>, rational<1, 2>), PowersList<Powers<1>, Powers<0>>{}));
    }

    SUBCASE("Test for a second order element")
    {
        constexpr auto form = PowersList<Powers<2>, Powers<1>, Powers<0>>{};

        constexpr auto constraints =
            std::tuple(EvaluateAt(rational<-1>), EvaluateAt(rational<0>), EvaluateAt(rational<1>));
        constexpr auto fns = derive_shape_functions(form, constraints);

        REQUIRE(
            get<0>(fns) == make_poly(
                               std::tuple(rational<0>, rational<-1, 2>, rational<1, 2>),
                               PowersList<Powers<0>, Powers<1>, Powers<2>>{}));

        REQUIRE(
            get<1>(fns) == make_poly(
                               std::tuple(rational<1>, rational<0>, rational<-1>),
                               PowersList<Powers<0>, Powers<1>, Powers<2>>{}));

        REQUIRE(
            get<2>(fns) == make_poly(
                               std::tuple(rational<0>, rational<1, 2>, rational<1, 2>),
                               PowersList<Powers<0>, Powers<1>, Powers<2>>{}));
    }

    SUBCASE("Test for a 3rd-order element with derivative DOF's")
    {
        constexpr auto form = PowersList<Powers<3>, Powers<2>, Powers<1>, Powers<0>>{};

        constexpr auto constraints = std::tuple(
            EvaluateAt(rational<-1>), EvaluateAt(rational<1>), Partial<0>::At(rational<-1>),
            Partial<0>::At(rational<1>));

        constexpr auto fns = derive_shape_functions(form, constraints);

        REQUIRE(
            fns[0] == make_poly(
                          std::tuple(rational<1, 2>, -rational<3, 4>, rational<0>, rational<1, 4>),
                          PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{}));

        REQUIRE(
            fns[1] == make_poly(
                          std::tuple(-rational<1, 4>, rational<0>, rational<3, 4>, rational<1, 2>),
                          PowersList<Powers<3>, Powers<2>, Powers<1>, Powers<0>>{}));

        REQUIRE(
            fns[2] == make_poly(
                          std::tuple(rational<1, 4>, -rational<1, 4>, -rational<1, 4>, rational<1, 4>),
                          PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{}));

        REQUIRE(
            fns[3] == make_poly(
                          std::tuple(-rational<1, 4>, -rational<1, 4>, rational<1, 4>, rational<1, 4>),
                          PowersList<Powers<0>, Powers<1>, Powers<2>, Powers<3>>{}));
    }
}

TEST_CASE("[Elements] Deriving bilinear shape functions")
{
    constexpr auto form = PowersList<Powers<1, 1>, Powers<1, 0>, Powers<0, 1>, Powers<0, 0>>{};
    constexpr auto constraints = std::tuple(
        EvaluateAt(rational<-1>, rational<-1>), EvaluateAt(rational<-1>, rational<1>),
        EvaluateAt(rational<1>, rational<1>), EvaluateAt(rational<1>, rational<-1>));

    constexpr auto fns = derive_shape_functions(form, constraints);

    REQUIRE(
        fns[0] == make_poly(
                      std::tuple(rational<1, 4>, -rational<1, 4>, -rational<1, 4>, rational<1, 4>),
                      PowersList<Powers<1, 1>, Powers<1, 0>, Powers<0, 1>, Powers<0, 0>>{}));

    REQUIRE(
        fns[1] == make_poly(
                      std::tuple(-rational<1, 4>, -rational<1, 4>, rational<1, 4>, rational<1, 4>),
                      PowersList<Powers<1, 1>, Powers<1, 0>, Powers<0, 1>, Powers<0, 0>>{}));

    REQUIRE(
        fns[2] == make_poly(
                      std::tuple(rational<1, 4>, rational<1, 4>, rational<1, 4>, rational<1, 4>),
                      PowersList<Powers<1, 1>, Powers<1, 0>, Powers<0, 1>, Powers<0, 0>>{}));
    
    REQUIRE(
        fns[3] == make_poly(
            std::tuple(-rational<1, 4>, rational<1, 4>, -rational<1, 4>, rational<1, 4>),
            PowersList<Powers<1, 1>, Powers<1, 0>, Powers<0, 1>, Powers<0, 0>>{}));
} // TEST_CASE

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 * End shape function API test.
 *******************************************************************************/
} // namespace Elements

} /* namespace Galerkin */

#endif /* ELEMENTS_HPP */

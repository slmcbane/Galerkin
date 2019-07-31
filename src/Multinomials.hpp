/*
 * Copyright (c) 2019, Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef MULTINOMIALS_HPP
#define MULTINOMIALS_HPP

#include "utils.hpp"
#include "Rationals.hpp"

namespace Galerkin
{

namespace Multinomials
{

using std::tuple;
using std::tuple_cat;

// I only support positive powers here.
template <unsigned v>
constexpr auto uint_constant = std::integral_constant<unsigned, v>();

// This struct represent the powers that each term is raised to in a multinomial.
template <unsigned... vs>
struct Powers
{
};

template <unsigned... vs>
constexpr auto nvars(Powers<vs...>) { return sizeof...(vs); }

template <auto I, unsigned... vs>
constexpr auto get_power(Powers<vs...>) { return std::get<I>(tuple(vs...)); }

template <unsigned... vs>
constexpr auto powers(std::integral_constant<unsigned, vs>...)
{
    return Powers<vs...>();
}

template <unsigned... vs, unsigned... ws>
constexpr bool operator<(Powers<vs...>, Powers<ws...>)
{
    static_assert(nvars(Powers<vs...>()) == nvars(Powers<ws...>()));

    return tuple(vs...) < tuple(ws...);
}

template <unsigned... vs, unsigned... ws>
constexpr bool operator==(Powers<vs...>, Powers<ws...>)
{
    static_assert(nvars(Powers<vs...>()) == nvars(Powers<ws...>()));

    return tuple(vs...) == tuple(ws...);
}

template <unsigned... vs, unsigned... ws>
constexpr bool operator<=(Powers<vs...>, Powers<ws...>)
{
    static_assert(nvars(Powers<vs...>()) == nvars(Powers<ws...>()));

    return tuple(vs...) <= tuple(ws...);
}

template <unsigned... vs, unsigned... ws>
constexpr auto operator+(Powers<vs...>, Powers<ws...>)
{
    static_assert(nvars(Powers<vs...>()) == nvars(Powers<ws...>()));
    return Powers<(ws + vs)...>();
}

// This struct represents a term in a multinomial.
// R is a Rational type and P is a Powers type.
template <class R, class P>
struct Term
{
};

template <class R, class P, class R2, class P2>
constexpr bool operator<(Term<R, P>, Term<R2, P2>)
{
    return P() < P2();
}

template <class R, class P, class R2, class P2>
constexpr bool operator==(Term<R, P>, Term<R2, P2>)
{
    return P() == P2();
}

template <class R, class P, class R2, class P2>
constexpr bool operator<=(Term<R, P>, Term<R2, P2>)
{
    return P() <= P2();
}

template <class R, unsigned... vs>
constexpr auto nvars(Term<R, Powers<vs...>>) { return nvars(Powers<vs...>()); }

template <int I, class R, unsigned... vs>
constexpr auto get_power(Term<R, Powers<vs...>>) { return get_power<I>(Powers<vs...>()); }

template <class R, unsigned... vs>
constexpr auto get_powers(Term<R, Powers<vs...>>) { return Powers<vs...>(); }

template <class R, class P>
constexpr auto coeff(Term<R, P>)
{
    return Rationals::rational<R::num(), R::den()>;
}

template <auto N, auto D, unsigned... vs>
constexpr auto term(Rationals::Rational<N, D>, Powers<vs...>)
{
    return Term<decltype(Rationals::rational<N, D>), Powers<vs...>>();
}

template <class R, class P>
constexpr auto operator-(Term<R, P>)
{
    return term(-R(), P());
}

// Helper to ensure that all terms in a multinomial have the same number of
// variables.
namespace
{
    constexpr auto map_nvars = [](auto Term) { return uint_constant<nvars(Term)>; };
}

// This struct is the type-level representation of a multinomial.
template <class... Terms>
struct Multinomial : public typeconst_list<Terms...>
{
    static_assert( sizeof...(Terms) == 0 ||
        typeconst_list<Terms...>().map(map_nvars).unique().count() == 1,
        "The number of variables in all terms of a Multinomial must be equal"
    );
};

template <class... T1s, class... T2s>
constexpr bool operator==(Multinomial<T1s...>, Multinomial<T2s...>)
{
    return std::is_same_v<Multinomial<T1s...>, Multinomial<T2s...>>;
}

template <auto I, class... Terms>
constexpr auto get_term(Multinomial<Terms...> mult)
{
    return get<I>(static_cast<typeconst_list<Terms...>>(mult));
}

template <class... Terms>
constexpr auto mult_from_list(typeconst_list<Terms...>)
{
    return Multinomial<Terms...>();
}

template <class... Terms>
constexpr auto combine_duplicate_powers(typeconst_list<Terms...> lst)
{
    if constexpr (sizeof...(Terms) == 1)
    {
        return mult_from_list(lst);
    }
    else
    {
        if constexpr (lst.head() == lst.tail().head())
        {
            constexpr auto c = coeff(lst.head()) + coeff(lst.tail().head());
            constexpr auto powers = get_powers(lst.head());
            return combine_duplicate_powers(
                make_list(term(c, powers)).append(lst.tail().tail()));
        }
        else
        {
            return mult_from_list(
                make_list(lst.head()).append(combine_duplicate_powers(lst.tail())));
        }
    }
}

template <>
constexpr auto combine_duplicate_powers(typeconst_list<>)
{
    return typeconst_list<>();
}

template <class... Terms>
constexpr auto drop_zeros(typeconst_list<Terms...> lst)
{
    if constexpr (coeff(lst.head()) == Rationals::rational<0>)
    {
        return drop_zeros(lst.tail());
    }
    else
    {
        return make_list(lst.head()).append(drop_zeros(lst.tail()));
    }
}

template <>
constexpr auto drop_zeros(typeconst_list<>)
{
    return typeconst_list<>();
}

template <class... Terms>
constexpr auto multinomial(Terms... ts)
{
    constexpr auto lst = Multinomial<Terms...>::sorted();
    return mult_from_list(drop_zeros(combine_duplicate_powers(lst)));
}

template <class... Terms>
constexpr auto nterms(Multinomial<Terms...>)
{
    return Multinomial<Terms...>::count();
}

/********************************************************************************
 * Test the basic constructor interface of a multinomial
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using Rationals::rational;

TEST_CASE("[Galerkin::Multinomials] Creating multinomials")
{
    SUBCASE("Test that powers get sorted")
    {
        constexpr auto mult = multinomial(
            term(rational<1>, powers(uint_constant<1>, uint_constant<1>)),
            term(rational<1>, powers(uint_constant<0>, uint_constant<1>)),
            term(rational<1>, powers(uint_constant<0>, uint_constant<0>)),
            term(rational<1>, powers(uint_constant<1>, uint_constant<0>)));

        REQUIRE(get_term<0>(mult) == term(rational<1>, powers(uint_constant<0>, uint_constant<0>)));
        REQUIRE(get_term<1>(mult) == term(rational<1>, powers(uint_constant<0>, uint_constant<1>)));
        REQUIRE(get_term<2>(mult) == term(rational<1>, powers(uint_constant<1>, uint_constant<0>)));
        REQUIRE(get_term<3>(mult) == term(rational<1>, powers(uint_constant<1>, uint_constant<1>)));
        REQUIRE(nterms(mult) == 4);
    }

    SUBCASE("Test that terms with the same power get combined")
    {
        constexpr auto mult = multinomial(
            term(rational<1, 2>, powers(uint_constant<0>)),
            term(rational<1>, powers(uint_constant<1>)),
            term(rational<1, 2>, powers(uint_constant<2>)),
            term(rational<1, 2>, powers(uint_constant<0>)));

        REQUIRE(nterms(mult) == 3);
        REQUIRE(get_term<0>(mult) == term(rational<1>, powers(uint_constant<0>)));
        REQUIRE(get_term<1>(mult) == term(rational<1>, powers(uint_constant<1>)));
        REQUIRE(get_term<2>(mult) == term(rational<1, 2>, powers(uint_constant<2>)));
    }

    SUBCASE("Test that terms with a zero coefficient get dropped")
    {
        constexpr auto mult = multinomial(
            term(rational<1>, powers(uint_constant<1>)),
            term(rational<1>, powers(uint_constant<0>)),
            term(rational<0>, powers(uint_constant<2>)));

        REQUIRE(nterms(mult) == 2);
        REQUIRE(get_term<0>(mult) == term(rational<1>, powers(uint_constant<0>)));
        REQUIRE(get_term<1>(mult) == term(rational<1>, powers(uint_constant<1>)));
    }
}

#endif

/********************************************************************************
 * End of basic construction of multinomial tests.
 *******************************************************************************/

template <class... Terms1, class... Terms2>
constexpr auto operator+(Multinomial<Terms1...> m1, Multinomial<Terms2...> m2)
{
    return multinomial(Terms1()..., Terms2()...);
}

template <class... Terms1, class... Terms2>
constexpr auto operator-(Multinomial<Terms1...> m1, Multinomial<Terms2...> m2)
{
    return multinomial(Terms1()..., -Terms2()...);
}

template <class... Terms>
constexpr auto operator-(Multinomial<Terms...>)
{
    return multinomial(-Terms()...);
}

/********************************************************************************
 * Test addition and subtraction of multinomials.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Multinomials] Addition and subtraction of multinomials")
{
    SUBCASE("All terms have the same powers")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1, 2>, powers(uint_constant<0>, uint_constant<0>)),
            term(rational<1, 3>, powers(uint_constant<0>, uint_constant<1>)),
            term(rational<1>, powers(uint_constant<1>, uint_constant<0>))
        );

        constexpr auto mult2 = multinomial(
            term(rational<1, 2>, powers(uint_constant<0>, uint_constant<0>)),
            term(rational<1, 3>, powers(uint_constant<0>, uint_constant<1>)),
            term(-rational<1, 2>, powers(uint_constant<1>, uint_constant<0>))
        );

        constexpr auto mult3 = multinomial(
            term(rational<1>, powers(uint_constant<0>, uint_constant<0>)),
            term(rational<2, 3>, powers(uint_constant<0>, uint_constant<1>)),
            term(rational<1, 2>, powers(uint_constant<1>, uint_constant<0>))
        );

        REQUIRE(mult1 + mult2 == mult3);
        REQUIRE(mult3 - mult2 == mult1);
        REQUIRE(mult3 - mult1 == mult2);
        REQUIRE(mult3 - mult1 - mult2 == multinomial());

        constexpr auto mult4 = multinomial(
            term(rational<3, 2>, powers(uint_constant<1>, uint_constant<0>))
        );

        REQUIRE(mult1 - mult2 == mult4);
        REQUIRE(mult4 + mult2 == mult1);
        REQUIRE(mult1 - mult2 - mult4 == multinomial());
    }

    SUBCASE("Terms have different powers, should merge")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(uint_constant<0>)),
            term(rational<1>, powers(uint_constant<2>))
        );

        constexpr auto mult2 = multinomial(
            term(rational<1>, powers(uint_constant<1>))
        );

        constexpr auto mult3 = multinomial(
            term(rational<1>, powers(uint_constant<0>)),
            term(rational<1>, powers(uint_constant<1>)),
            term(rational<1>, powers(uint_constant<2>))
        );

        constexpr auto mult4 = multinomial(
            term(rational<1>, powers(uint_constant<0>)),
            term(rational<1>, powers(uint_constant<2>)),
            term(-rational<1>, powers(uint_constant<1>))
        );

        REQUIRE(mult1 + mult2 == mult3);
        REQUIRE(mult3 - mult2 == mult1);
        REQUIRE(mult3 - mult1 == mult2);
        REQUIRE(-mult1 + -mult2 == -mult3);

        REQUIRE(mult1 - mult2 == mult4);
    }
}

#endif

/********************************************************************************
 * End of multinomial addition and subtraction tests.
 *******************************************************************************/

template <class R1, class P1, class R2, class P2>
constexpr auto operator*(Term<R1, P1>, Term<R2, P2>)
{
    return term(R1()*R2(), P1()+P2());
}

template <class R, class P, class... Terms>
constexpr auto operator*(Term<R, P> t, Multinomial<Terms...>)
{
    return multinomial(t * Terms()...);
}

template <class... Terms1, class... Terms2>
constexpr auto operator*(Multinomial<Terms1...> mult1, Multinomial<Terms2...> mult2)
{
    return static_sum<0, mult2.count()>(
        [=](auto I) { return get_term<I()>(mult2) * mult1; },
        multinomial()
    );
}

/********************************************************************************
 * Test multiplication of multinomials (by multinomial and by scalar constant).
 * Division by a scalar is implemented, but not multinomial division.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Multinomials] Test multiplication and division of multinomials")
{
    SUBCASE("Multiplication of polynomials")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(uint_constant<2>)),
            term(rational<1, 2>, powers(uint_constant<0>))
        );

        constexpr auto mult2 = multinomial(
            term(rational<1>, powers(uint_constant<2>)),
            term(-rational<1, 2>, powers(uint_constant<0>))
        );

        constexpr auto mult3 = multinomial(
            term(rational<1>, powers(uint_constant<4>)),
            term(-rational<1, 4>, powers(uint_constant<0>))
        );

        REQUIRE(mult1 * mult2 == mult3);
    }
}

#endif

/********************************************************************************
 * End multiplication and division tests
 *******************************************************************************/

} /* namespace Multinomials */

} /* namespace Galerkin */

#endif /* MULTINOMIALS_HPP */

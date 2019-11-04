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

// This struct represent the powers that each term is raised to in a multinomial.
template <auto... vs>
struct Powers
{
};

template <class... Args, auto... vs>
constexpr auto raise(std::tuple<Args...> args, Powers<vs...>)
{
    static_assert(sizeof...(Args) == sizeof...(vs));
    static_assert(sizeof...(Args) > 0);
    return static_reduce<1, sizeof...(Args), 1>(
        [=]([[maybe_unused]] auto I) {
            constexpr auto power = std::get<I()>(std::tuple(vs...));
            auto arg = std::get<I()>(args);
            return static_reduce<0, power, 1>(
                [=]([[maybe_unused]] auto I) { return arg; },
                Rationals::rational<1>,
                [] (auto x, auto y) { return x * y; }
            );
        },
        static_reduce<0, std::get<0>(std::tuple(vs...)), 1>(
            [=]([[maybe_unused]] auto I) { return std::get<0>(args); },
            Rationals::rational<1>,
            [] (auto x, auto y) { return x * y; }
        ),
        [] (auto x, auto y) { return x * y; }
    );
}

template <auto... vs>
constexpr auto nvars(Powers<vs...>) { return sizeof...(vs); }

template <auto I, auto... vs>
constexpr auto get_power(Powers<vs...>) { return std::get<I>(std::tuple(vs...)); }

template <auto... vs>
constexpr auto powers(std::integral_constant<decltype(vs), vs>...)
{
    return Powers<vs...>();
}

template <auto... vs, auto... ws>
constexpr bool operator<(Powers<vs...>, Powers<ws...>)
{
    static_assert(nvars(Powers<vs...>()) == nvars(Powers<ws...>()));

    return std::tuple(vs...) < std::tuple(ws...);
}

template <auto... vs, auto... ws>
constexpr bool operator==(Powers<vs...>, Powers<ws...>)
{
    static_assert(nvars(Powers<vs...>()) == nvars(Powers<ws...>()));

    return std::tuple(vs...) == std::tuple(ws...);
}

template <auto v, auto... vs, auto w, auto... ws>
constexpr bool operator<=(Powers<v, vs...>, Powers<w, ws...>)
{
    static_assert(nvars(Powers<v, vs...>()) == nvars(Powers<w, ws...>()));

    if constexpr (sizeof...(vs) == 0)
    {
        return v <= w;
    }
    else if constexpr (v == w)
    {
        return Powers<vs...>() <= Powers<ws...>();
    }
    else if constexpr (v < w)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <auto... vs, auto... ws>
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
    template <class... Args>
    constexpr auto operator()(std::tuple<Args...> args)
    {
        return R() * raise(args, P());
    }
};

template <class R, class P, class R2, class P2>
constexpr bool operator<(Term<R, P>, Term<R2, P2>)
{
    return P() < P2();
}

template <class R, class P, class R2, class P2>
constexpr bool operator==(Term<R, P>, Term<R2, P2>)
{
    return std::is_same_v<P, P2> && std::is_same_v<R, R2>;
}

template <class R, class P, class R2, class P2>
constexpr bool operator<=(Term<R, P>, Term<R2, P2>)
{
    return P() <= P2();
}

template <class R, auto... vs>
constexpr auto nvars(Term<R, Powers<vs...>>) { return nvars(Powers<vs...>()); }

template <int I, class R, auto... vs>
constexpr auto get_power(Term<R, Powers<vs...>>) { return get_power<I>(Powers<vs...>()); }

template <class R, auto... vs>
constexpr auto get_powers(Term<R, Powers<vs...>>) { return Powers<vs...>(); }

template <class R, class P>
constexpr auto coeff(Term<R, P>)
{
    return Rationals::rational<R::num(), R::den()>;
}

template <auto N, auto D, auto... vs>
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
constexpr auto map_nvars = [](auto Term) { return intgr_constant<nvars(Term)>; };
}

// This struct is the type-level representation of a multinomial.
template <class... Terms>
struct Multinomial : public typeconst_list<Terms...>
{
    static_assert(sizeof...(Terms) == 0 ||
                      typeconst_list<Terms...>().map(map_nvars).unique().count() == 1,
                  "The number of variables in all terms of a Multinomial must be equal");

    template <class... Types>
    constexpr auto operator()(Types... args) const noexcept
    {
        auto tup = std::tuple(args...);
        return this->operator()(tup);
    }

    template <class... Types>
    constexpr auto operator()(std::tuple<Types...> args) const noexcept
    {
        if constexpr (sizeof...(Terms) == 0)
        {
            return intgr_constant<0>;
        }
        else
        {
            return static_sum<1, sizeof...(Terms)>(
                [=](auto I) {
                    return get<I()>(static_cast<typeconst_list<Terms...>>(*this))(args);
                },
                this->head()(args));
        }
    }
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
        if constexpr (get_powers(lst.head()) == get_powers(lst.tail().head()))
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
constexpr auto multinomial(Terms...)
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
            term(rational<1>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<0>, intgr_constant<0>)),
            term(rational<1>, powers(intgr_constant<1>, intgr_constant<0>)));

        REQUIRE(get_term<0>(mult) == term(rational<1>, powers(intgr_constant<0>, intgr_constant<0>)));
        REQUIRE(get_term<1>(mult) == term(rational<1>, powers(intgr_constant<0>, intgr_constant<1>)));
        REQUIRE(get_term<2>(mult) == term(rational<1>, powers(intgr_constant<1>, intgr_constant<0>)));
        REQUIRE(get_term<3>(mult) == term(rational<1>, powers(intgr_constant<1>, intgr_constant<1>)));
        REQUIRE(nterms(mult) == 4);
    }

    SUBCASE("Test that terms with the same power get combined")
    {
        constexpr auto mult = multinomial(
            term(rational<1, 2>, powers(intgr_constant<0>)),
            term(rational<1>, powers(intgr_constant<1>)),
            term(rational<1, 2>, powers(intgr_constant<2>)),
            term(rational<1, 2>, powers(intgr_constant<0>)));

        REQUIRE(nterms(mult) == 3);
        REQUIRE(get_term<0>(mult) == term(rational<1>, powers(intgr_constant<0>)));
        REQUIRE(get_term<1>(mult) == term(rational<1>, powers(intgr_constant<1>)));
        REQUIRE(get_term<2>(mult) == term(rational<1, 2>, powers(intgr_constant<2>)));
    }

    SUBCASE("Test that terms with a zero coefficient get dropped")
    {
        constexpr auto mult = multinomial(
            term(rational<1>, powers(intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<0>)),
            term(rational<0>, powers(intgr_constant<2>)));

        REQUIRE(nterms(mult) == 2);
        REQUIRE(get_term<0>(mult) == term(rational<1>, powers(intgr_constant<0>)));
        REQUIRE(get_term<1>(mult) == term(rational<1>, powers(intgr_constant<1>)));
    }
}

#endif

/********************************************************************************
 * End of basic construction of multinomial tests.
 *******************************************************************************/

template <class... Terms1, class... Terms2>
constexpr auto operator+(Multinomial<Terms1...>, Multinomial<Terms2...>)
{
    return multinomial(Terms1()..., Terms2()...);
}

template <class... Terms1, class... Terms2>
constexpr auto operator-(Multinomial<Terms1...>, Multinomial<Terms2...>)
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
            term(rational<1, 2>, powers(intgr_constant<0>, intgr_constant<0>)),
            term(rational<1, 3>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<1>, intgr_constant<0>)));

        constexpr auto mult2 = multinomial(
            term(rational<1, 2>, powers(intgr_constant<0>, intgr_constant<0>)),
            term(rational<1, 3>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(-rational<1, 2>, powers(intgr_constant<1>, intgr_constant<0>)));

        constexpr auto mult3 = multinomial(
            term(rational<1>, powers(intgr_constant<0>, intgr_constant<0>)),
            term(rational<2, 3>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1, 2>, powers(intgr_constant<1>, intgr_constant<0>)));

        REQUIRE(mult1 + mult2 == mult3);
        REQUIRE(mult3 - mult2 == mult1);
        REQUIRE(mult3 - mult1 == mult2);
        REQUIRE(mult3 - mult1 - mult2 == multinomial());

        constexpr auto mult4 = multinomial(
            term(rational<3, 2>, powers(intgr_constant<1>, intgr_constant<0>)));

        REQUIRE(mult1 - mult2 == mult4);
        REQUIRE(mult4 + mult2 == mult1);
        REQUIRE(mult1 - mult2 - mult4 == multinomial());
    }

    SUBCASE("Terms have different powers, should merge")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(intgr_constant<0>)),
            term(rational<1>, powers(intgr_constant<2>)));

        constexpr auto mult2 = multinomial(
            term(rational<1>, powers(intgr_constant<1>)));

        constexpr auto mult3 = multinomial(
            term(rational<1>, powers(intgr_constant<0>)),
            term(rational<1>, powers(intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<2>)));

        constexpr auto mult4 = multinomial(
            term(rational<1>, powers(intgr_constant<0>)),
            term(rational<1>, powers(intgr_constant<2>)),
            term(-rational<1>, powers(intgr_constant<1>)));

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
    return term(R1() * R2(), P1() + P2());
}

template <class R1, class P1, class I, I v>
constexpr auto operator*(Term<R1, P1>, std::integral_constant<I, v>)
{
    return term(R1() * Rationals::rational<v>, P1());
}

template <class R1, class P1, class I, I v>
constexpr auto operator/(Term<R1, P1>, std::integral_constant<I, v>)
{
    return term(R1() / Rationals::rational<v>, P1());
}

template <class R1, class P1, class I, I v>
constexpr auto operator*(std::integral_constant<I, v>, Term<R1, P1>)
{
    return term(R1() * Rationals::rational<v>, P1());
}

template <class R, class P, class... Terms>
constexpr auto operator*(Term<R, P> t, Multinomial<Terms...>)
{
    return multinomial(t * Terms()...);
}

template <class R, class P, auto N, auto D>
constexpr auto operator*(Term<R, P>, Rationals::Rational<N, D>)
{
    return term(Rationals::rational<N, D> * R(), P());
}

template <class R, class P, auto N, auto D>
constexpr auto operator/(Term<R, P>, Rationals::Rational<N, D>)
{
    return term(R() / Rationals::rational<N, D>, P());
}

template <class R, class P, auto N, auto D>
constexpr auto operator*(Rationals::Rational<N, D>, Term<R, P>)
{
    return term(Rationals::rational<N, D> * R(), P());
}

template <class... Terms1, class... Terms2>
constexpr auto operator*(Multinomial<Terms1...> mult1, Multinomial<Terms2...> mult2)
{
    return static_sum<0, mult2.count()>(
        [=](auto I) { return get_term<I()>(mult2) * mult1; },
        multinomial());
}

template <class... Terms, class I, I v>
constexpr auto operator*(Multinomial<Terms...>, std::integral_constant<I, v>)
{
    return multinomial(std::integral_constant<I, v>() * Terms()...);
}

template <class... Terms, auto N, auto D>
constexpr auto operator*(Multinomial<Terms...>, Rationals::Rational<N, D>)
{
    return multinomial(Rationals::rational<N, D> * Terms()...);
}

template <class... Terms, class I, I v>
constexpr auto operator/(Multinomial<Terms...>, std::integral_constant<I, v>)
{
    return multinomial( (Terms() / std::integral_constant<I, v>())... );
}

template <class... Terms, auto N, auto D>
constexpr auto operator/(Multinomial<Terms...>, Rationals::Rational<N, D>)
{
    return multinomial( (Terms() / Rationals::rational<N, D>)... );
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
            term(rational<1>, powers(intgr_constant<2>)),
            term(rational<1, 2>, powers(intgr_constant<0>)));

        constexpr auto mult2 = multinomial(
            term(rational<1>, powers(intgr_constant<2>)),
            term(-rational<1, 2>, powers(intgr_constant<0>)));

        constexpr auto mult3 = multinomial(
            term(rational<1>, powers(intgr_constant<4>)),
            term(-rational<1, 4>, powers(intgr_constant<0>)));

        REQUIRE(mult1 * mult2 == mult3);

        constexpr auto mult4 = multinomial(
            term(rational<1>, powers(intgr_constant<2>)),
            term(rational<2>, powers(intgr_constant<1>)),
            term(rational<3>, powers(intgr_constant<0>)));

        constexpr auto mult5 = multinomial(
            term(rational<1>, powers(intgr_constant<4>)),
            term(rational<2>, powers(intgr_constant<3>)),
            term(rational<7, 2>, powers(intgr_constant<2>)),
            term(rational<1>, powers(intgr_constant<1>)),
            term(rational<3, 2>, powers(intgr_constant<0>)));

        REQUIRE(mult1 * mult4 == mult5);
    }

    SUBCASE("Multiplication of multinomials")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(intgr_constant<2>, intgr_constant<0>)),
            term(rational<2>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<3>, powers(intgr_constant<0>, intgr_constant<0>)));

        constexpr auto mult2 = multinomial(
            term(rational<3>, powers(intgr_constant<0>, intgr_constant<2>)),
            term(rational<2>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<0>, intgr_constant<0>)));

        constexpr auto mult3 = multinomial(
            term(rational<3>, powers(intgr_constant<2>, intgr_constant<2>)),
            term(rational<2>, powers(intgr_constant<2>, intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<2>, intgr_constant<0>)),
            term(rational<6>, powers(intgr_constant<1>, intgr_constant<2>)),
            term(rational<4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(rational<2>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<9>, powers(intgr_constant<0>, intgr_constant<2>)),
            term(rational<6>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<3>, powers(intgr_constant<0>, intgr_constant<0>)));

        REQUIRE(mult1 * mult2 == mult3);
    }

    SUBCASE("Multiplication by a scalar")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(intgr_constant<4>)),
            term(rational<2>, powers(intgr_constant<3>)),
            term(rational<7, 2>, powers(intgr_constant<2>)),
            term(rational<1>, powers(intgr_constant<1>)),
            term(rational<3, 2>, powers(intgr_constant<0>)));

        constexpr auto mult2 = multinomial(
            term(rational<1, 2>, powers(intgr_constant<4>)),
            term(rational<1>, powers(intgr_constant<3>)),
            term(rational<7, 4>, powers(intgr_constant<2>)),
            term(rational<1, 2>, powers(intgr_constant<1>)),
            term(rational<3, 4>, powers(intgr_constant<0>))
        );

        constexpr auto mult3 = multinomial(
            term(rational<-3>, powers(intgr_constant<4>)),
            term(rational<-6>, powers(intgr_constant<3>)),
            term(-rational<21, 2>, powers(intgr_constant<2>)),
            term(-rational<3>, powers(intgr_constant<1>)),
            term(-rational<9, 2>, powers(intgr_constant<0>))
        );

        REQUIRE(mult1 * rational<1, 2> == mult2);
        REQUIRE(mult1 * std::integral_constant<int, -3>() == mult3);
    }

    SUBCASE("Division by a scalar")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(intgr_constant<2>, intgr_constant<0>)),
            term(rational<2>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<3>, powers(intgr_constant<0>, intgr_constant<0>)));

        REQUIRE(mult1 / std::integral_constant<int, -3>() == mult1 * -rational<1, 3>);
        REQUIRE(mult1 / rational<3, 2> == mult1 * rational<2, 3>);
    }
}

#endif

/********************************************************************************
 * End multiplication and division tests
 *******************************************************************************/

/********************************************************************************
 * Test application of a multinomial
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Multinomials] Multinomial application")
{
    SUBCASE("Applying a polynomial")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(intgr_constant<2>)),
            term(rational<2>, powers(intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<0>)));

        REQUIRE(mult1(rational<1, 2>) == rational<9, 4>);
        REQUIRE(mult1(0.5) == doctest::Approx(2.25));
        REQUIRE(mult1(0.5f) == doctest::Approx(2.25f));
    }

    SUBCASE("Applying a true multinomial")
    {
        constexpr auto mult = multinomial(
            term(rational<3, 4>, powers(intgr_constant<2>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(rational<1, 2>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<1, 3>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<2, 3>, powers(intgr_constant<0>, intgr_constant<0>))
        );

        REQUIRE(mult(rational<0>, rational<0>) == rational<2, 3>);
        REQUIRE(mult(rational<1>, rational<1>) == rational<5, 2>);
        REQUIRE(mult(rational<1, 2>, rational<1, 3>) == rational<163, 144>);
        REQUIRE(mult(0.5, 1.0 / 3) == doctest::Approx(163.0 / 144));
    }
}

#endif

/********************************************************************************
 * End of multinomial application tests.
 *******************************************************************************/

template <auto... vs>
constexpr auto powers_from_tuple(std::tuple<std::integral_constant<decltype(vs), vs>...>)
{
    return Powers<vs...>();
}

template <auto I, auto v, auto... vs>
constexpr auto subtract_one(
    std::tuple<std::integral_constant<decltype(v), v>,
    std::integral_constant<decltype(vs), vs>...>)
{
    if constexpr (I == 0)
    {
        return std::tuple_cat(
            std::tuple(intgr_constant<v-1>),
            std::make_tuple(intgr_constant<vs>...)
        );
    }
    else
    {
        return std::tuple_cat(
            std::tuple(intgr_constant<v>),
            subtract_one<I-1>(std::make_tuple(intgr_constant<vs>...))
        );
    }
}

template <auto I, auto... vs>
constexpr auto subtract_one(Powers<vs...>)
{
    return powers_from_tuple(subtract_one<I>(std::tuple(intgr_constant<vs>...)));
}

template <auto I, class R, class P>
constexpr auto partial(Term<R, P>)
{
    static_assert(I < nvars(P()));
    constexpr auto pow = get_power<I>(P());
    if constexpr (pow == 0)
    {
        return term(Rationals::rational<0>, P());
    }
    else
    {
        return term(R() * intgr_constant<get_power<I>(P())>, subtract_one<I>(P()));
    }
}

template <auto I, class... Terms>
constexpr auto partial(Multinomial<Terms...>)
{
    static_assert(I < nvars(get_term<0>(Multinomial<Terms...>())));
    return multinomial(partial<I>(Terms())...);
}

/********************************************************************************
 * Test for partial derivatives of a multinomial
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Partial derivatives of a multinomial")
{
    SUBCASE("Test for a polynomial")
    {
        constexpr auto mult1 = multinomial(
            term(rational<1>, powers(intgr_constant<2>)),
            term(rational<2>, powers(intgr_constant<1>)),
            term(rational<1>, powers(intgr_constant<0>)));

        constexpr auto mult2 = multinomial(
            term(rational<2>, powers(intgr_constant<1>)),
            term(rational<2>, powers(intgr_constant<0>))
        );

        REQUIRE(partial<0>(mult1) == mult2);
    }

    SUBCASE("Test for a multinomial")
    {
        constexpr auto mult = multinomial(
            term(rational<3, 4>, powers(intgr_constant<2>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(rational<1, 2>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<1, 3>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<2, 3>, powers(intgr_constant<0>, intgr_constant<0>))
        );

        constexpr auto dmult0 = multinomial(
            term(rational<6, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1, 2>, powers(intgr_constant<0>, intgr_constant<0>))
        );

        constexpr auto dmult1 = multinomial(
            term(rational<3, 4>, powers(intgr_constant<2>, intgr_constant<0>)),
            term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<1, 3>, powers(intgr_constant<0>, intgr_constant<0>))
        );

        REQUIRE(partial<0>(mult) == dmult0);
        REQUIRE(partial<1>(mult) == dmult1);
    }
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 * End of partial derivative tests.
 *******************************************************************************/

} /* namespace Multinomials */

} /* namespace Galerkin */

#endif /* MULTINOMIALS_HPP */

#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include "Multinomials.hpp"

namespace Galerkin
{

namespace Legendre
{

template <auto I>
constexpr auto polynomial = 
    (Multinomials::multinomial(
        term(Rationals::rational<2*I-1>,
             Multinomials::powers(Multinomials::uint_constant<1>))) * polynomial<I-1>
    - polynomial<I-2> * Multinomials::uint_constant<I-1>) /
    Multinomials::uint_constant<I>;

template <>
constexpr auto polynomial<0> = Multinomials::multinomial(
    term(Rationals::rational<1>,
         Multinomials::powers(Multinomials::uint_constant<0>))
);

template <>
constexpr auto polynomial<1> = Multinomials::multinomial(
    term(Rationals::rational<1>,
         Multinomials::powers(Multinomials::uint_constant<1>))
);

/********************************************************************************
 * Test that Legendre polynomials are correct.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using namespace Multinomials;
using namespace Rationals;

TEST_CASE("Test computed Legendre polynomials")
{
    REQUIRE(polynomial<0> == Multinomials::multinomial(term(rational<1>, powers(uint_constant<0>))));
    REQUIRE(polynomial<1> == Multinomials::multinomial(term(rational<1>, powers(uint_constant<1>))));
    REQUIRE(polynomial<2> ==
        Multinomials::multinomial(
            term(rational<3, 2>, powers(uint_constant<2>)),
            term(-rational<1, 2>, powers(uint_constant<0>))
        )
    );

    REQUIRE(polynomial<10> ==
        multinomial(
            term(rational<46189, 256>, powers(uint_constant<10>)),
            term(-rational<109395, 256>, powers(uint_constant<8>)),
            term(rational<90090, 256>, powers(uint_constant<6>)),
            term(-rational<30030, 256>, powers(uint_constant<4>)),
            term(rational<3465, 256>, powers(uint_constant<2>)),
            term(-rational<63, 256>, powers(uint_constant<0>))
        )
    );
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

} /* namespace Legendre */

} /* namespace Galerkin */

#endif /* LEGENDRE_HPP */
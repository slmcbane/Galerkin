#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include "ElementBase.hpp"
#include "Elements.hpp"

namespace Galerkin
{

namespace Elements
{

template <int Degree, class T>
class IntervalElement : public ElementBase<IntervalElement<Degree, T>>
{
    T m_scaling;
    T m_translation;
public:
    constexpr static auto basis = derive_shape_functions(
        powers_up_to(intgr_constant<Degree>),
        evenly_spaced(Rationals::rational<-1>, Rationals::rational<1>, intgr_constant<Degree>)
            .map([](auto N) { return EvaluateAt<decltype(N)>{}; })
    );

    IntervalElement(T a, T b) : m_scaling((b-a) / 2), m_translation((a+b) / 2)
    {}
};

} // namespace Elements

template <int Degree, class T>
struct DefaultIntegrationOrder<Elements::IntervalElement<Degree, T>>
{
    constexpr static int order = 2 * Degree;
};

namespace Elements
{

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::Elements] Check IntervalElement basis functions")
{
    IntervalElement<1, double> elt(0.0, 1.0);

    REQUIRE(get<0>(elt.basis) ==
        Metanomials::metanomial(Metanomials::term(Rationals::rational<-1, 2>,
            Metanomials::Powers<1>{}), Metanomials::term(Rationals::rational<1, 2>,
            Metanomials::Powers<0>{})));
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace Elements

} /* namespace Galerkin */

#endif /* INTERVAL_HPP */
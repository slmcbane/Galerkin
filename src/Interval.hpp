#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include "ElementBase.hpp"
#include "Elements.hpp"
#include "UniformScaling.hpp"

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

    constexpr IntervalElement(T a, T b) : m_scaling((b-a) / 2), m_translation((a+b) / 2)
    {}

    constexpr auto coordinate_map() const noexcept
    {
        return Transforms::UniformScaling(m_scaling, std::array<T, 1>{m_translation});
    }
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

    IntervalElement<2, double> elt2(0.0, 1.0);

    REQUIRE(get<1>(elt2.basis) ==
        Metanomials::metanomial(
            Metanomials::term(Rationals::rational<-1>, Metanomials::Powers<2>{}),
            Metanomials::term(Rationals::rational<1>, Metanomials::Powers<0>{})
        ));
} // TEST_CASE

TEST_CASE("[Galerkin::Elements] Check mapping of points through IntervalElement")
{
    IntervalElement<1, double> elt(0.5, 1.0);

    REQUIRE(elt.transform(std::tuple(0.0))[0] == doctest::Approx(0.75));
    REQUIRE(elt.transform(std::array<float, 1>{0.5})[0] == doctest::Approx(0.875));
} // TEST_CASE

TEST_CASE("[Galerkin::Elements] Test computed mass and stiffness matrices for IntervalElement")
{

SUBCASE("No symmetric tag")
{
    constexpr IntervalElement<2, double> elt(0.5, 1.0);

    constexpr auto mass_matrix = elt.form_matrix(
        [](auto f, auto g) { return f * g; }
    );

    REQUIRE(mass_matrix(0, 0) == doctest::Approx(1.0 / 15));
    REQUIRE(mass_matrix(0, 1) == doctest::Approx(1.0 / 30));
    REQUIRE(mass_matrix(0, 2) == doctest::Approx(-1.0 / 60));
    REQUIRE(mass_matrix(1, 1) == doctest::Approx(4.0 / 15));
    REQUIRE(mass_matrix(1, 2) == doctest::Approx(1.0 / 30));
    REQUIRE(mass_matrix(2, 2) == doctest::Approx(1.0 / 15));
    REQUIRE(mass_matrix(0, 1) == doctest::Approx(mass_matrix(1, 0)));
    REQUIRE(mass_matrix(0, 2) == doctest::Approx(mass_matrix(2, 0)));
    REQUIRE(mass_matrix(1, 2) == doctest::Approx(mass_matrix(2, 1)));
} // SUBCASE

} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

} // namespace Elements

} /* namespace Galerkin */

#endif /* INTERVAL_HPP */
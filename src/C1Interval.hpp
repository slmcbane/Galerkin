#ifndef GALERKIN_C1INTERVAL_HPP
#define GALERKIN_C1INTERVAL_HPP

#include "ElementBase.hpp"
#include "Elements.hpp"
#include "Polynomial.hpp"
#include "UniformScaling.hpp"

namespace Galerkin
{

namespace Elements
{

template <class T>
class C1IntervalElement : public ElementBase<C1IntervalElement<T>>
{
    T m_scaling;
    T m_translation;

  public:
    constexpr static auto base_basis = derive_shape_functions(
        Polynomials::PowersList<
            Polynomials::Powers<0>, Polynomials::Powers<1>, Polynomials::Powers<2>,
            Polynomials::Powers<3>>{},
        std::tuple(
            Elements::Partial<0>::At<Rationals::Rational>(Rationals::rational<-1>),
            Elements::EvaluateAt(Rationals::rational<-1>), Elements::EvaluateAt(Rationals::rational<1>),
            Elements::Partial<0>::At<Rationals::Rational>(Rationals::rational<1>)));

    constexpr static auto basis_size = base_basis.size();

    constexpr auto basis() const noexcept
    {
        using result_type = decltype(base_basis[0] * m_scaling);
        return std::array<result_type, 4>{
            m_scaling * base_basis[0], one<T> * base_basis[1], one<T> * base_basis[2],
            m_scaling * base_basis[3]};
    }

    constexpr auto coordinate_map() const noexcept
    {
        return Transforms::UniformScaling(m_scaling, std::array<T, 1>{m_translation});
    }

    constexpr C1IntervalElement(const T &s, const T &t) noexcept : m_scaling(s), m_translation(t) {}
};

} // namespace Elements

template <class T>
struct DefaultIntegrationOrder<Elements::C1IntervalElement<T>>
{
    constexpr static int order = 6;
};

} // namespace Galerkin

#endif // GALERKIN_C1INTERVAL_HPP

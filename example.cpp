#include "src/ElementBase.hpp"
#include "src/Elements.hpp"
#include "src/TriangleTransform.hpp"

using Galerkin::Elements::evaluate_at;
using Galerkin::Metanomials::Powers;
using Galerkin::Rationals::rational;

class TetElement;

// This template specialization must come *before* the definition of the
// TetElement class because the instantiation of ElementBase depends on the
// specialization of DefaultIntegrationOrder.
template <>
class Galerkin::DefaultIntegrationOrder<TetElement> : public Galerkin::IntegrationOrder<2>
{
};

class TetElement : public Galerkin::Elements::ElementBase<TetElement>
{
public:
    template <class... Points>
    constexpr TetElement(const Points &... points) noexcept : m_trans(points...)
    {
    }

    static constexpr auto basis = Galerkin::Elements::derive_shape_functions(
        Galerkin::Elements::make_form(
            Powers<1, 0>{}, Powers<0, 1>{}, Powers<0, 0>{}),
        Galerkin::make_list(
            evaluate_at(rational<-1>, rational<-1>),
            evaluate_at(rational<-1>, rational<1>),
            evaluate_at(rational<1>, rational<-1>)));

    constexpr auto &coordinate_map() const noexcept { return m_trans; }

private:
    Galerkin::Transforms::TriangleTransform<double> m_trans;
};

extern "C"
{
    void mass_matrix(const double *verts, double *dst);
    void stiffness_matrix(const double *verts, double *dst);
} // extern "C"

void mass_matrix(const double *verts, double *dst)
{
    const auto p1 = std::tuple(verts[0], verts[1]);
    const auto p2 = std::tuple(verts[2], verts[3]);
    const auto p3 = std::tuple(verts[4], verts[5]);
    const TetElement element(p1, p2, p3);

    constexpr auto form = [](auto f, auto g) { return f * g; };

    const auto matrix = element.form_matrix(form);
    int index = 0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i; j < 3; ++j)
        {
            dst[index++] = matrix(i, j);
        }
    }
}

void stiffness_matrix(const double *verts, double *dst)
{
    const auto p1 = std::tuple(verts[0], verts[1]);
    const auto p2 = std::tuple(verts[2], verts[3]);
    const auto p3 = std::tuple(verts[4], verts[5]);
    const TetElement element(p1, p2, p3);

    auto form = [&](auto f, auto g) {
        return element.template partial<0>(f) * element.template partial<0>(g) +
               element.template partial<1>(f) * element.template partial<1>(g);
    };

    // Gradient of the first order shape functions is constant, so zero-th order
    // integration is fine.
    const auto matrix = element.form_matrix(form, Galerkin::IntegrationOrder<0>{});
    int index = 0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i; j < 3; ++j)
        {
            dst[index++] = matrix(i, j);
        }
    }
}